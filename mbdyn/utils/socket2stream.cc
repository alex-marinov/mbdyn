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
#include <iostream>

#include "sock.h"

void
usage(int rc)
{
	silent_cerr(
"MBDyn (C) is a multibody analysis code.\n"
"http://www.mbdyn.org\n"
"\n"
"Copyright (C) 1996-2004\n"
"\n"
"Pierangelo Masarati	<masarati@aero.polimi.it>\n"
"\n"
"Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano\n"
"via La Masa, 34 - 20156 Milano, Italy\n"
"http://www.aero.polimi.it\n"
"\n"
"usage: socket2stream [options]\n"
"\n"
"    -h <host>\t"	"host name (for INET sockets; default: \"localhost\")\n"
"    -n <channels>\t"	"number of channels (default: auto-detect)\n"
"    -p <port>\t"	"port number (for INET sockets)\n"
"    -P <path>\t"	"path (for LOCAL sockets)\n"
"    -s\t\t"		"decrease verbosity level\n"
			<< std::endl);
	exit(rc);
}

int
main(int argc, char *argv[])
{
	int	nChannels = 0;
	char	*host = 0;
	int	port = -1;
	char	*path = 0;
	char	*next;

	while (true) {
		int opt = getopt(argc, argv, "h:n:p:P:s");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'h':
			host = optarg;
			break;

		case 'n':
			nChannels = strtol(optarg, &next, 10);
			if (next[0] != '\0') {
				silent_cerr("unable to parse option -n "
						"\"" << optarg << "\""
						<< std::endl);
				usage(EXIT_FAILURE);
			}
			break;

		case 'p':
			port = strtol(optarg, &next, 10);
			if (next[0] != '\0') {
				silent_cerr("unable to parse option -p "
						"\"" << optarg << "\""
						<< std::endl);
				usage(EXIT_FAILURE);
			}
			break;

		case 'P':
			path = optarg;
			break;

		case 's':
			::fSilent++;
			break;

		default:
			usage(EXIT_SUCCESS);
		}
	}

	if (path == 0 && (host == 0 && port == -1)) {
		usage(EXIT_FAILURE);
	}
}

