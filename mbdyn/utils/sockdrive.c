/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

#include <stdio.h>
#include <stdlib.h>

#ifdef USE_SOCKET_DRIVES

#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#include <string.h>
#include <signal.h>

#include "sock.h"

const unsigned short int PORT = 5555;
const char *SERVERHOST = "localhost";
const char *SERVERPATH = "/var/mbdyn/mbdyn.sock";

static void
usage(void)
{
	fprintf(stderr,
		"\n\tusage: sockdrive [h:p:D:w:Wi:I:] label [value]\n\n"
		"\t\t-h host\t\thost name\n"
		"\t\t-p port\t\tport number\n"
		"\t\t-P path\t\tpath of named pipe\n"
		"\t\t-D user\t\tuser name\n"
		"\t\t-w cred\t\tuser credentials\n"
		"\t\t-W\t\tprompt for user credentials\n"
		"\t\t-i {yes|no}\tincremental (i.e. value[label] += value)\n"
		"\t\t-I {yes|no}\timpulsive (reset after one step)\n\n"
		"\tlabel:\tfile drive index to modify\n"
		"\tvalue:\tnew value (or increment if -i)\n\n");
}

int
main(int argc, char *argv[])
{
	int sock;
	struct sockaddr_in server_name;

	char *path = NULL;
	char *host = (char *)SERVERHOST;
	unsigned short int port = PORT;

	char *user = NULL;
	char *cred = NULL;

	int inc = 0;
	int imp = 0;
	char *label;
	char *value = NULL;

	FILE *fd;


	while (1) {
		int opt;

		opt = getopt (argc, argv, "h:p:P:D:w:Wi:I:");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'h':
			host = optarg;
			break;

		case 'p':
			port = atoi(optarg);
			break;

		case 'P':
			path = optarg;
			break;

		case 'D':
			user = optarg;
			break;

		case 'w':
			cred = strdup(optarg);
			break;

		case 'W':
			if (cred) {
				free(cred);
			}

			if (user) {
				char buf[1024];

				snprintf(buf, sizeof(buf), "Password for user \"%s\": ", user);
				cred = getpass(buf);

			} else {
				cred = getpass("Password: ");
			}

			if (cred) {
				cred = strdup(cred);
			}
			break;

		case 'i':
			if (strcasecmp(optarg, "yes") == 0) {
				inc = 1;
			} else if (strcasecmp(optarg, "no") == 0) {
				inc = -1;
			}
			break;

		case 'I':
			if (strcasecmp(optarg, "yes") == 0) {
				imp = 1;
			} else if (strcasecmp(optarg, "no") == 0) {
				imp = -1;
			}
			break;
		}
	}

	if ((argc - optind) < 1) {
		usage();
		exit(EXIT_SUCCESS);
	}
	label = argv[optind];

	if ((argc - optind) > 1) {
		value = argv[optind + 1];
	}

	if (path) {
		sock = make_named_socket(path, 0);
	} else {
		sock = make_inet_socket(&server_name, host, port, 0);
	}
	if (sock < 0) {
		fprintf(stderr, "socket initialization error\n");
		exit(EXIT_FAILURE);
	}

	fd = fdopen(sock, "w");
	if (user) {
		fprintf(fd, "user: %s\n", user);
		if (cred) {
			fprintf(fd, "password: %s\n", cred);
			memset(cred, '\0', strlen(cred));
		}
	}

	fprintf(fd, "label: %s\n", label);

	if (inc == 1) {
		fprintf(fd, "inc: yes\n");
	} else if (inc == -1) {
		fprintf(fd, "inc: no\n");
	}

	if (imp == 1) {
		fprintf(fd, "imp: yes\n");
	} else if (imp == -1) {
		fprintf(fd, "imp: no\n");
	}

	if (value != NULL) {
		fprintf(fd, "value: %s\n", value);
	}

	fprintf(fd, ".\n");

	fclose(fd);

	exit(EXIT_SUCCESS);
}

#else /* ! USE_SOCKET_DRIVES */

int
main(void)
{
	fprintf(stderr, "need sockets\n");
	exit(EXIT_FAILURE);
}

#endif /* ! USE_SOCKET_DRIVES */

