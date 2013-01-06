/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#ifdef USE_SOCKET

#include <cstring>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include "ac/getopt.h"

#include <iostream>
#include <cstring>
#include <sstream>
#include <vector>

#include "sock.h"
#include "s2s.h"

static struct s2s_t	*s2s_list = 0;

s2s_t::s2s_t(void)
: nChannels(0),
path(0),
host(0),
port(-1),
create(false),
block(BLOCK_UNKNOWN),
sock(-1),
stream2socket(true),
progname(0)
{
	buf[0] = '\0';

	next = ::s2s_list;
	::s2s_list = next;
}

s2s_t::~s2s_t(void)
{
	struct s2s_t	**p;

	for (p = &::s2s_list; *p != 0; p = &(*p)->next) {
		if (*p == this) {
			*p = (*p)->next;
			break;
		}
	}
}

void
s2s_t::usage(int rc) const
{
	const char
			*b = "",
			*C = "",
			*F = "",
			*h = "",
			*n = "",
			*p = "",
			*P = "",
			*s = "";

	b = "    -b {y|n}\t\t"		"blocking mode (default: yes)\n";
	C = "    -C\t\t\t"		"create socket (default: connect to existing)\n";
	h = "    -h <host>\t\t"		"host name (for INET sockets; default: \"localhost\")\n";
	n = "    -n <channels>\t"	"number of channels (default: auto-detect)\n";
	p = "    -p <port>\t\t"		"port number (for INET sockets)\n";
	P = "    -P <path>\t\t"		"path (for LOCAL sockets)\n";
	s = "    -s\t\t\t"		"decrease verbosity level\n";

	if (!this->stream2socket) {
		F = "    -F <format>\t\t"	"output format (%[.<precision>]{eEfF})\n";
	}

	silent_cerr(
"\n"
"    MBDyn (C) is a multibody analysis code.\n"
"    http://www.mbdyn.org\n"
"\n"
"    Copyright (C) 1996-2013\n"
"\n"
"    Pierangelo Masarati	<masarati@aero.polimi.it>\n"
"\n"
"    Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano\n"
"    via La Masa, 34 - 20156 Milano, Italy\n"
"    http://www.aero.polimi.it\n"
"\n"
"    usage: " << this->progname << " [options]\n"
"\n"
			<< b
			<< C
			<< F
			<< h
			<< n
			<< p
			<< P
			<< s
			<< std::endl);
	exit(rc);
}

void
s2s_t::shutdown(void)
{
	if (this->sock >= 0) {
		::shutdown(this->sock, SHUT_RDWR);
		sock = -1;
	}

	if (this->path) {
		unlink(this->path);
		this->path = 0;
	}
}

void
s2s_shutdown(int signum)
{
	struct s2s_t	**p, **nextp;

	for (p = &::s2s_list; *p != 0; p = nextp ) {
		nextp = &(*p)->next;
		(*p)->shutdown();
		*p = 0;
	}

	::s2s_list = 0;

	exit(EXIT_SUCCESS);
}

static void
parse_format(const char *fmt)
{
	char		*f = (char *)fmt, *next;

	if (f[0] == '%') {
		f++;
	}

	if (f[0] == '.') {
		f++;
		unsigned precision = strtoul(f, &next, 10);
		if (next == f) {
			silent_cerr("unable to parse \"precision\" "
				"in format string "
				"\"" << fmt << "\"" << std::endl);
			throw;
		}
		std::cout.precision(precision);
		f = next;
	}
	
	if (f[0] == 'e' || f[0] == 'E') {
		std::cout.setf(std::ios::scientific);

	} else if (f[0] != 'f' && f[0] != 'F') {
		silent_cerr("unable to parse format string "
			"\"" << fmt << "\"" << std::endl);
		throw;
	}
}

void
s2s_t::parse(int argc, char *argv[])
{
	char	*next;

	next = std::strrchr(argv[0], '/');
	if (next != 0) {
		next++;

	} else {
		next = argv[0];
	}

	// libtool paranoia
	if (strncmp(next, "lt-", STRLENOF("lt-")) == 0) {
		next += STRLENOF("lt-");
	}

	if (strcmp(next, "socket2stream") == 0) {
		this->stream2socket = false;
		this->progname = "socket2stream";

	} else if (strcmp(next, "stream2socket") == 0) {
		this->stream2socket = true;
		this->progname = "stream2socket";

	} else {
		silent_cerr("tool name mismatch: \"" << argv[0] << "\"" << std::endl);
		exit(EXIT_FAILURE);
	}

	const char	*optstring;
	if (this->stream2socket) {
		optstring = "b:Ch:n:p:P:s";

	} else {
		optstring = "b:CF:h:n:p:P:s";
	}

	while (true) {
		int opt = getopt(argc, argv, optstring);

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'b':
			if (strcmp(optarg, "y") == 0) {
				this->block = BLOCK_YES;
			} else if (strcmp(optarg, "n") == 0) {
				this->block = BLOCK_NO;
			} else {
				silent_cerr("invalid value "
						"\"" << optarg << "\""
						" for option -b"
						<< std::endl);
				usage(EXIT_FAILURE);
			}
			break;

		case 'C':
			this->create = true;
			break;

		case 'F':
			parse_format(optarg);
			break;

		case 'h':
			this->host = optarg;
			break;

		case 'n':
			this->nChannels = strtoul(optarg, &next, 10);
			if (next[0] != '\0') {
				silent_cerr("unable to parse option -n "
						"\"" << optarg << "\""
						<< std::endl);
				usage(EXIT_FAILURE);
			}
			break;

		case 'p':
			this->port = strtoul(optarg, &next, 10);
			if (next[0] != '\0') {
				silent_cerr("unable to parse option -p "
						"\"" << optarg << "\""
						<< std::endl);
				usage(EXIT_FAILURE);
			}
			break;

		case 'P':
			this->path = optarg;
			break;

		case 's':
			::fSilent++;
			break;

		default:
			usage(EXIT_SUCCESS);
		}
	}

	if (this->path == 0 && (this->host == 0 && this->port == -1)) {
		usage(EXIT_FAILURE);
	}
}

void
s2s_t::prepare(void)
{
	s2s_sockaddr_t	addr;
	int		save_errno;
	struct sockaddr	*addrp = 0;
	socklen_t	addrl = 0;

	if (this->path) {
		addr.ms_type = AF_LOCAL;
		addr.ms_domain = PF_LOCAL;
		addrp = (struct sockaddr *)&addr.ms_addr.ms_addr_local;
		addrl = sizeof(addr.ms_addr.ms_addr_local);
		strncpy(this->buf, this->path, sizeof(this->buf));
		
		if (this->create) {
			this->sock = mbdyn_make_named_socket(&addr.ms_addr.ms_addr_local,
				this->path, 1, &save_errno);
		
			if (this->sock == -1) {
				const char	*err_msg = strerror(save_errno);

      				silent_cerr("socket(" << this->buf << ") failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
      				throw;

 		  	} else if (this->sock == -2) {
				const char	*err_msg = strerror(save_errno);

		      		silent_cerr("bind(" << this->buf << ") failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
		      		throw;
		   	}

		} else {
			addr.ms_addr.ms_addr_local.sun_family = AF_UNIX;
			strncpy(addr.ms_addr.ms_addr_local.sun_path, this->path,
					sizeof(addr.ms_addr.ms_addr_local.sun_path));
		}

	} else {
		addr.ms_type = AF_INET;
		addr.ms_domain = PF_INET;
		addrp = (struct sockaddr *)&addr.ms_addr.ms_addr_inet;
		addrl = sizeof(addr.ms_addr.ms_addr_inet);

		if (this->host == 0) {
			this->host = "127.0.0.1";

		} else {
			char *p = std::strchr(const_cast<char *>(this->host), ':');
			if (p != 0) {
				char	*next;

				p[0] = '\0';
				p++;
				this->port = strtoul(p, &next, 10);
				if (next[0] != '\0') {
					silent_cerr("unable to parse port out of "
						"\"" << this->host << ":" << p << "\""
						<< std::endl);
					exit(EXIT_FAILURE);
				}
			}
		}
		snprintf(this->buf, sizeof(this->buf), "%s:%u", this->host, this->port);

		if (this->create) {
	 	  	this->sock = mbdyn_make_inet_socket(&addr.ms_addr.ms_addr_inet, 
					this->host, this->port, 1, &save_errno);

			if (this->sock == -1) {
				const char	*err_msg = strerror(save_errno);

      				silent_cerr("socket(" << this->buf << ") failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
				throw;

   			} else if (this->sock == -2) {
				const char	*err_msg = strerror(save_errno);

      				silent_cerr("bind(" << this->buf << ") failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
				throw;

			} else if (this->sock == -3) {
     				silent_cerr("illegal host[:port] name \"" << this->buf << "\" "
					"(" << save_errno << ")"
					<< std::endl);
				throw;
  			}

		} else {
			addr.ms_addr.ms_addr_inet.sin_family = AF_INET;
			addr.ms_addr.ms_addr_inet.sin_port = htons(this->port);
					
			if (inet_aton(this->host, &addr.ms_addr.ms_addr_inet.sin_addr) == 0) {
				silent_cerr("unknown host \"" << this->host << "\"" << std::endl);
				throw;	
			}
		}
	}

	if (this->create) {
		if (listen(this->sock, 1) < 0) {
			save_errno = errno;
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("listen(" << this->sock << "," << this->buf << ") failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw;
   		}

		int sock = this->sock;
		this->sock = accept(sock, addrp, &addrl);
		close(sock);
		if (this->sock == -1) {
			save_errno = errno;
			const char	*err_msg = strerror(save_errno);

                	silent_cerr("accept(" << this->sock << ",\"" << this->buf << "\") "
				"failed (" << save_errno << ": " << err_msg << ")"
				<< std::endl);
			throw;
        	}

	} else {
		this->sock = socket(addr.ms_domain, SOCK_STREAM, 0);
		if (this->sock < 0){
			save_errno = errno;
			const char	*err_msg = strerror(save_errno);
				
			silent_cerr("socket(" << this->buf << ") failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
			throw;
		}

		if (connect(this->sock, addrp, addrl) < 0) {
			save_errno = errno;
			const char	*err_msg = strerror(save_errno);
				
			silent_cerr("connect(" << this->sock << ",\"" << this->buf << "\"," << addrl << ") "
				"failed (" << save_errno << ": " << err_msg << ")"
				<< std::endl);
			throw;
		}
	}

	signal(SIGTERM, s2s_shutdown);
	signal(SIGINT, s2s_shutdown);
	signal(SIGPIPE, s2s_shutdown);
}

bool
s2s_t::is_blocking(void) const
{
	return block == BLOCK_YES;
}

ssize_t
s2s_t::send(int flags) const
{
	switch (block) {
	case BLOCK_NO:
		flags |= MSG_DONTWAIT;
		break;

	case BLOCK_YES:
		flags &= ~MSG_DONTWAIT;
		break;
	}

	return ::send(sock, (char *)&dbuf[0], sizeof(double)*nChannels, flags);
}

ssize_t
s2s_t::recv(int flags)
{
	switch (block) {
	case BLOCK_NO:
		flags |= MSG_DONTWAIT;
		break;

	case BLOCK_YES:
		flags &= ~MSG_DONTWAIT;
		break;
	}

	return ::recv(sock, (char *)&dbuf[0], sizeof(double)*nChannels, flags);
}

#endif // USE_SOCKET
