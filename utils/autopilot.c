/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include <termios.h>
#include "ac/getopt.h"

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <sys/un.h>

#include <string.h>

#ifdef HAVE_SASL2
#if defined(HAVE_SASL_SASL_H)
#include <sasl/sasl.h>
#elif defined (HAVE_SASL_H)
#include <sasl.h>
#endif /* HAVE_SASL_SASL_H || HAVE_SASL_H */
#include "mbsasl.h"
#endif /* HAVE_SASL2 */

#include "sock.h"

const unsigned short int	PORT = 5555;
const char			*SERVERHOST = "localhost";
const char			*SERVERPATH = "/var/mbdyn/mbdyn.sock";

static void
keys(FILE * fh)
{
	fprintf(fh,
			"\tkeys:\n"
			"\t\t'i':\tswitches the incremental mode on\n"
			"\t\t'p':\tincrements the drive\n"
			"\t\t'm':\tdecrements the drive\n"
			"\t\t'^D':\tquits\n\n");
}

static void
usage(void)
{
 	fprintf(stderr,
			"\n\tusage: autopilot [h:p:D:vw:Wx:] <label>\n\n"
			"\t\t-D <user>\tuser name\n"
			"\t\t-h <host>\thost name\n"
			"\t\t-m <mech>\tSASL mechanism(s)\n"
			"\t\t-p <port>\tport number\n"
			"\t\t-P <path>\tpath for named socked\n"
			"\t\t-S\t\tuse SASL"
#ifndef HAVE_SASL2
				" (not supported)"
#endif /* ! HAVE_SASL2 */
			"\n"
			"\t\t-v\t\tverbose\n"
			"\t\t-w <cred>\tuser credentials\n"
			"\t\t-W\t\tprompt for user credentials\n"
			"\t\t-x <value>\tincrement\n"
			"\n"
			"\t\t<label>:\tfile drive (base 1) index to modify\n\n");
	keys(stderr);
}

/* Use this variable to remember original terminal attributes. */
struct termios saved_attributes;

void
reset_input_mode(void)
{
	tcsetattr(STDIN_FILENO, TCSANOW, &saved_attributes);
}

void
set_input_mode(void)
{
	struct termios tattr;

	/* Make sure stdin is a terminal. */
	if (!isatty(STDIN_FILENO)) {
		fprintf(stderr, "Not a terminal.\n");
		exit(EXIT_FAILURE);
	}

	/* Save the terminal attributes so we can restore them later. */
	tcgetattr(STDIN_FILENO, &saved_attributes);
	atexit(reset_input_mode);

	/* Set the funny terminal modes. */
	tcgetattr(STDIN_FILENO, &tattr);
	tattr.c_lflag &= ~(ICANON|ECHO); /* Clear ICANON and ECHO. */
	tattr.c_cc[VMIN] = 0;
	tattr.c_cc[VTIME] = 0;
	if (tcsetattr(STDIN_FILENO, TCSAFLUSH, &tattr) < 0) {
		perror("tcsetattr");
		exit(EXIT_FAILURE);
	}
}

static int sasl = 0;
#ifdef HAVE_SASL2
static struct mbdyn_sasl_t mbdyn_sasl = MBDYN_SASL_INIT;
#endif /* HAVE_SASL2 */

const char		*path = NULL;
const char		*host = NULL;
unsigned short int	port = 0;

int
send_message(const char *message)
{
	int sock;
	struct sockaddr_in peer_name;
	FILE *fd;

	/* Create the socket. */
	if (path) {
		sock = mbdyn_make_named_socket(0, path, 0, NULL);
	} else {
		sock = mbdyn_make_inet_socket(0, host, port, 0, NULL);
	}
	if (sock < 0) {
		return -1;
	}

	/* Connect to the server. */
	if (0 > connect(sock, (struct sockaddr *) &peer_name,
				sizeof(peer_name))) {
		return -1;
	}

	if (sasl) {
#ifdef HAVE_SASL2
		if (mbdyn_sasl_client_auth(sock, NULL, &mbdyn_sasl) != SASL_OK) {
			return -1;
		}
#endif /* HAVE_SASL2 */
	}

	fd = fdopen(sock, "w");
	fputs(message, fd);
	fclose(fd);

	return 0;
}

int
main(int argc, char *argv[])
{
	char *user = NULL;
	char *cred = NULL;

	char *increment = "1.";
	char *label = NULL;

	char c;

	char *auth = NULL;
	char *inc = NULL;
	char *plus = NULL;
	char *minus = NULL;

	char *mech = NULL;

	int verbose = 0;

	host = SERVERHOST;
	port = PORT;

	while (1) {
		int opt;

		opt = getopt(argc, argv, "D:h:m:p:P:Svw:Wx:");

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

		case 'v':
			verbose++;
			break;

		case 'w':
			cred = strdup(optarg);
			break;

		case 'W': {
			char *tmp = getpass("password: ");

			if (tmp) {
				cred = strdup(tmp);
				memset(tmp, '\0', strlen(tmp));
			}
			break;
		}

		case 'x':
			increment = optarg;
			break;
		}
	}

	if (argc - optind < 1) {
		usage();
		exit(EXIT_SUCCESS);
	}

	label = argv[optind];

	/* messaggi: */
	if (!sasl) {
		if (user) {
			size_t l = strlen(user);
			if (cred) {
				l += strlen(cred) + 7 + 11;

				auth = (char *)calloc(sizeof(char), l + 1);
				snprintf(auth, l + 1, "user: %s\npassword: %s\n", user, cred);
			} else {
				l += 7;

				auth = (char *)calloc(sizeof(char), l + 1);
				snprintf(auth, l + 1, "user: %s\n", user);
			}
		}
	} else {
#ifdef HAVE_SASL2
		if (verbose) {
			printf("initializing SASL data...\n");
		}

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
#endif /* HAVE_SASL2 */
	}

	if (auth) {
		size_t l = strlen(auth) + 8 + strlen(label) + 9 + 2;

		inc = (char *)calloc(sizeof(char), l + 1);
		snprintf(inc, l, "%slabel: %s\ninc: yes\n.\n", auth, label);

		l = strlen(auth) + 8 + strlen(label) + 8 + strlen(increment) + 2;
		plus = (char *)calloc(sizeof(char), l + 1);
		minus = (char *)calloc(sizeof(char), l + 1 + 1);
		snprintf(plus, l + 1, "%slabel: %s\nvalue: %s\n.\n", auth, label, increment);
		snprintf(minus, l + 2, "%slabel: %s\nvalue: -%s\n.\n", auth, label, increment);

	} else {
		size_t l = 8 + strlen(label) + 9 + 2;
		
		inc = (char *)calloc(sizeof(char), l + 1);
		snprintf(inc, l + 1, "label: %s\ninc: yes\n.\n", label);

		l = 8 + strlen(label) + 8 + strlen(increment) + 2;
		plus = (char *)calloc(sizeof(char), l + 1);
		minus = (char *)calloc(sizeof(char), l + 1 + 1);
		snprintf(plus, l + 1, "label: %s\nvalue: %s\n.\n", label, increment);
		snprintf(minus, l + 2, "label: %s\nvalue: -%s\n.\n", label, increment);
	}

	set_input_mode();

	if (verbose) {
		fprintf(stdout, "Connecting to host %s:%d\n", host, port);
		if (user) {
			if (cred) {
				fprintf(stdout, "Accounting as \"%s\" (with creds)\n", user);
			
			} else {
				fprintf(stdout, "Accounting as \"%s\"\n", user);
			}
		}
		fprintf(stdout, "Incrementing drive %s by %s\n", label, increment);
		keys(stdout);
	}

	while (1) {
		size_t i;

		i = read(STDIN_FILENO, &c, 1);

		if (i > 0) {
			char	*msg = NULL;

			if (c == '\004') {         /* `C-d' */
				break;
			}

			switch (c) {
			case 'i':
				msg = inc;
				break;

				break;

			case 'p':
				msg = plus;
				break;

			case 'm':
				msg = minus;
				break;

			}

			if (send_message(msg) == -1) {
				fprintf(stderr, "unable to connect to host %s:%d\n",
						host, port);
			}
		}
	}

	if (cred) {
		free(cred);
	}

	exit(EXIT_SUCCESS);
}

