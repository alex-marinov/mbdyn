/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <stdio.h>
#include <stdlib.h>
#include "ac/getopt.h"

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h>

#include <string.h>
#include <signal.h>

#ifdef HAVE_SASL2
#if defined(HAVE_SASL_SASL_H)
#include <sasl/sasl.h>
#elif defined (HAVE_SASL_H)
#include <sasl.h>
#endif /* HAVE_SASL_SASL_H || HAVE_SASL_H */
#include "mbsasl.h"
#endif /* HAVE_SASL2 */

#include "sock.h"

const unsigned short int PORT = 5555;
const char *SERVERHOST = "localhost";
const char *SERVERPATH = "/var/mbdyn/mbdyn.sock";

static void
usage(void)
{
	fprintf(stderr,
		"\n\tusage: sockdrive [h:p:D:w:Wi:I:] <label> [<value>]\n\n"
		"\t\t-D <user>\tuser name\n"
		"\t\t-h <host>\thost name\n"
		"\t\t-i {yes|no}\tincremental (i.e. value[<label>] += <value>)\n"
		"\t\t-I {yes|no}\timpulsive (reset after one step)\n"
		"\t\t-m <mech>\tSASL mechanism (needs -S)\n"
		"\t\t-p <port>\tport number\n"
		"\t\t-P <path>\tpath of named pipe\n"
		"\t\t-S\t\tenable SASL auth"
#ifndef HAVE_SASL2
			" (not supported)"
#endif /* ! HAVE_SASL2 */
		"\n"
		"\t\t-w <cred>\tuser credentials\n"
		"\t\t-W\t\tprompt for user credentials\n"
		"\t\t<label>:\tfile drive (base 1) index to modify\n"
		"\n"
		"\t\t<value>:\tnew value (or increment if -i)\n\n");
}

	static int sasl = 0;
#ifdef HAVE_SASL2
static struct mbdyn_sasl_t mbdyn_sasl = MBDYN_SASL_INIT;
#endif /* HAVE_SASL2 */

int
main(int argc, char *argv[])
{
	int sock;

	char *path = NULL;
	char *host = (char *)SERVERHOST;
	unsigned short int port = PORT;

	char *user = NULL;
	char *cred = NULL;
	char *mech = NULL;

	int inc = 0;
	int imp = 0;
	char *label;
	char *value = NULL;

	FILE *fd;


	while (1) {
		int opt;

		opt = getopt (argc, argv, "D:h:i:I:m:p:P:S:w:W");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'D':
			user = optarg;
			break;

		case 'h':
			host = optarg;
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

		case 'm':
			mech = optarg;
#ifndef HAVE_SASL2
			fprintf(stderr, "SASL not supported\n");
#endif /* ! HAVE_SASL2 */
			break;

		case 'p':
			port = atoi(optarg);
			break;

		case 'P':
			path = optarg;
			break;

		case 'S':
			sasl++;
#ifndef HAVE_SASL2
			fprintf(stderr, "SASL not supported\n");
			exit(EXIT_FAILURE);
#endif /* ! HAVE_SASL2 */
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
		sock = mbdyn_make_named_socket(0, path, 0, NULL);
	} else {
		sock = mbdyn_make_inet_socket(0, host, port, 0, NULL);
	}
	if (sock < 0) {
		fprintf(stderr, "socket initialization error\n");
		exit(EXIT_FAILURE);
	}

	if (sasl) {
#ifdef HAVE_SASL2
		printf("initializing SASL data...\n");

		mbdyn_sasl.use_sasl = MBDYN_SASL_CLIENT;
		mbdyn_sasl.sasl_flags = MBDYN_SASL_FLAG_CRITICAL | MBDYN_SASL_FLAG_USERAUTHZ;
		mbdyn_sasl.sasl_mech = mech;
		mbdyn_sasl.sasl_user = user;
		mbdyn_sasl.sasl_cred = cred;
		mbdyn_sasl.sasl_hostname = host;

		if (mbdyn_sasl_client_init(&mbdyn_sasl) != SASL_OK) {
			fprintf(stderr, "SASL init failed\n");
			exit(EXIT_FAILURE);
		}

		if (mbdyn_sasl_client_auth(sock, NULL, &mbdyn_sasl) != SASL_OK) {
			fprintf(stderr, "SASL auth failed\n");
			exit(EXIT_FAILURE);
		}
#endif /* HAVE_SASL2 */
	}

	fd = fdopen(sock, "w");

	if (!sasl) {
		if (user) {
			fprintf(fd, "user: %s\n", user);
			if (cred) {
				fprintf(fd, "password: %s\n", cred);
				memset(cred, '\0', strlen(cred));
			}
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

