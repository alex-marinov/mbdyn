/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <ac/getopt.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "sock.h"
#include "s2s.h"

s2s_t::s2s_t(void)
: sock(-1),
nChannels(0),
path(0),
host(0),
port(-1),
create(false),
stream2socket(true),
progname(0)
{
	buf[0] = '\0';
}

s2s_t	s2s;

static void
usage(int rc)
{
	char	*C = "",
		*F = "",
		*h = "",
		*n = "",
		*p = "",
		*P = "",
		*s = "";
	
	C = "    -C\t\t\t"		"create socket (default: connect to existing)\n";
	h = "    -h <host>\t\t"		"host name (for INET sockets; default: \"localhost\")\n";
	n = "    -n <channels>\t"	"number of channels (default: auto-detect)\n";
	p = "    -p <port>\t\t"		"port number (for INET sockets)\n";
	P = "    -P <path>\t\t"		"path (for LOCAL sockets)\n";
	s = "    -s\t\t\t"		"decrease verbosity level\n";

	if (!::s2s.stream2socket) {
		F = "    -F <format>\t\t"	"output format (%[.<precision>]{eEfF})\n";
	}

	silent_cerr(
"\n"
"    MBDyn (C) is a multibody analysis code.\n"
"    http://www.mbdyn.org\n"
"\n"
"    Copyright (C) 1996-2005\n"
"\n"
"    Pierangelo Masarati	<masarati@aero.polimi.it>\n"
"\n"
"    Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano\n"
"    via La Masa, 34 - 20156 Milano, Italy\n"
"    http://www.aero.polimi.it\n"
"\n"
"    usage: " << s2s.progname << " [options]\n"
"\n"
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
s2s_shutdown(int signum)
{
	if (::s2s.sock >= 0) {
		shutdown(::s2s.sock, SHUT_RDWR);
		::s2s.sock = -1;
	}

	if (::s2s.path) {
		unlink(::s2s.path);
		::s2s.path = 0;
	}

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

	} else if (f[0] != 'f' & f[0] != 'F') {
		silent_cerr("unable to parse format string "
			"\"" << fmt << "\"" << std::endl);
		throw;
	}
}

void
s2s_prepare(int argc, char *argv[])
{
	char	*next;

	next = strrchr(argv[0], '/');
	if (next != 0) {
		next++;

	} else {
		next = argv[0];
	}

	if (strcmp(next, "socket2stream") == 0) {
		::s2s.stream2socket = false;
		::s2s.progname = "socket2stream";

	} else if (strcmp(next, "stream2socket") == 0) {
		::s2s.stream2socket = true;
		::s2s.progname = "stream2socket";

	} else {
		silent_cerr("tool name mismatch: \"" << argv[0] << "\"" << std::endl);
		exit(EXIT_FAILURE);
	}

	char	*optstring;
	if (::s2s.stream2socket) {
		optstring = "Ch:n:p:P:s";

	} else {
		optstring = "CF:h:n:p:P:s";
	}

	while (true) {
		int opt = getopt(argc, argv, optstring);

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'C':
			::s2s.create = true;
			break;

		case 'F':
			parse_format(optarg);
			break;

		case 'h':
			::s2s.host = optarg;
			break;

		case 'n':
			::s2s.nChannels = strtoul(optarg, &next, 10);
			if (next[0] != '\0') {
				silent_cerr("unable to parse option -n "
						"\"" << optarg << "\""
						<< std::endl);
				usage(EXIT_FAILURE);
			}
			break;

		case 'p':
			::s2s.port = strtoul(optarg, &next, 10);
			if (next[0] != '\0') {
				silent_cerr("unable to parse option -p "
						"\"" << optarg << "\""
						<< std::endl);
				usage(EXIT_FAILURE);
			}
			break;

		case 'P':
			::s2s.path = optarg;
			break;

		case 's':
			::fSilent++;
			break;

		default:
			usage(EXIT_SUCCESS);
		}
	}

	if (::s2s.path == 0 && (::s2s.host == 0 && ::s2s.port == -1)) {
		usage(EXIT_FAILURE);
	}

	s2s_sockaddr_t	addr;
	int		save_errno;
	struct sockaddr	*addrp = 0;
	socklen_t	addrl = 0;

	if (::s2s.path) {
		addr.ms_type = AF_LOCAL;
		addr.ms_domain = PF_LOCAL;
		addrp = (struct sockaddr *)&addr.ms_addr.ms_addr_local;
		addrl = sizeof(addr.ms_addr.ms_addr_local);
		strncpy(s2s.buf, ::s2s.path, sizeof(s2s.buf));
		
		if (::s2s.create) {
			::s2s.sock = mbdyn_make_named_socket(::s2s.path, 1, &save_errno);
		
			if (::s2s.sock == -1) {
				const char	*err_msg = strerror(save_errno);

      				silent_cerr("socket(" << s2s.buf << ") failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
      				throw;

 		  	} else if (::s2s.sock == -2) {
				const char	*err_msg = strerror(save_errno);

		      		silent_cerr("bind(" << s2s.buf << ") failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
		      		throw;
		   	}

		} else {
			addr.ms_addr.ms_addr_local.sun_family = AF_UNIX;
			strncpy(addr.ms_addr.ms_addr_local.sun_path, ::s2s.path,
					sizeof(addr.ms_addr.ms_addr_local.sun_path));
		}

	} else {
		addr.ms_type = AF_INET;
		addr.ms_domain = PF_INET;
		addrp = (struct sockaddr *)&addr.ms_addr.ms_addr_inet;
		addrl = sizeof(addr.ms_addr.ms_addr_inet);

		if (::s2s.host == 0) {
			::s2s.host = "127.0.0.1";

		} else {
			char	*p = strchr(::s2s.host, ':');
			if (p != 0) {
				p[0] = '\0';
				p++;
				::s2s.port = strtoul(p, &next, 10);
				if (next[0] != '\0') {
					silent_cerr("unable to parse port out of "
						"\"" << ::s2s.host << ":" << p << "\""
						<< std::endl);
					exit(EXIT_FAILURE);
				}
			}
		}
		snprintf(s2s.buf, sizeof(s2s.buf), "%s:%u", ::s2s.host, ::s2s.port);

		if (::s2s.create) {
	 	  	::s2s.sock = mbdyn_make_inet_socket(&addr.ms_addr.ms_addr_inet, 
					::s2s.host, ::s2s.port, 1, &save_errno);

			if (::s2s.sock == -1) {
				const char	*err_msg = strerror(save_errno);

      				silent_cerr("socket(" << s2s.buf << ") failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
				throw;

   			} else if (::s2s.sock == -2) {
				const char	*err_msg = strerror(save_errno);

      				silent_cerr("bind(" << s2s.buf << ") failed "
					"(" << save_errno << ": " << err_msg << ")"
					<< std::endl);
				throw;

			} else if (::s2s.sock == -3) {
     				silent_cerr("illegal host[:port] name \"" << s2s.buf << "\" "
					"(" << save_errno << ")"
					<< std::endl);
				throw;
  			}

		} else {
			addr.ms_addr.ms_addr_inet.sin_family = AF_INET;
			addr.ms_addr.ms_addr_inet.sin_port = htons(::s2s.port);
					
			if (inet_aton(::s2s.host, &addr.ms_addr.ms_addr_inet.sin_addr) == 0) {
				silent_cerr("unknown host \"" << ::s2s.host << "\"" << std::endl);
				throw;	
			}
		}
	}

	if (::s2s.create) {
		if (listen(::s2s.sock, 1) < 0) {
			save_errno = errno;
			const char	*err_msg = strerror(save_errno);

      			silent_cerr("listen(" << ::s2s.sock << "," << s2s.buf << ") failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
      			throw;
   		}

		int sock = ::s2s.sock;
		::s2s.sock = accept(sock, addrp, &addrl);
		close(sock);
		if (::s2s.sock == -1) {
			save_errno = errno;
			const char	*err_msg = strerror(save_errno);

                	silent_cerr("accept(" << ::s2s.sock << ",\"" << s2s.buf << "\") "
				"failed (" << save_errno << ": " << err_msg << ")"
				<< std::endl);
			throw;
        	}

	} else {
		::s2s.sock = socket(addr.ms_domain, SOCK_STREAM, 0);
		if (::s2s.sock < 0){
			save_errno = errno;
			const char	*err_msg = strerror(save_errno);
				
			silent_cerr("socket(" << s2s.buf << ") failed "
				"(" << save_errno << ": " << err_msg << ")"
				<< std::endl);
			throw;
		}

		if (connect(::s2s.sock, addrp, addrl) < 0) {
			save_errno = errno;
			const char	*err_msg = strerror(save_errno);
				
			silent_cerr("connect(" << ::s2s.sock << ",\"" << s2s.buf << "\"," << addrl << ") "
				"failed (" << save_errno << ": " << err_msg << ")"
				<< std::endl);
			throw;
		}
	}

	signal(SIGTERM, s2s_shutdown);
	signal(SIGINT, s2s_shutdown);
	signal(SIGPIPE, s2s_shutdown);
}

