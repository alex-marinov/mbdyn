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
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include "ac/getopt.h"

#include <string.h>

#include "crypt.h"

static void
usage(int rc)
{
	fprintf(stderr, "usage: crypt [-f salt] [-h] [-c asserted] cred\n");
	exit(rc);
}

int
main(int argc, char *argv[])
{
#ifndef HAVE_CRYPT
	fprintf(stderr, "no valid crypt available\n");
	exit(EXIT_FAILURE);
#else /* HAVE_CRYPT */

	char	*salt_format = NULL, salt[35];
	char	*asserted_cred = NULL, asserted_cred_buf[35];
	size_t	asserted_cred_len = 0;

	while (1) {
		int	opt = getopt(argc, argv, "c:f:h");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'c':
			asserted_cred = optarg;
			break;

		case 'f':
			salt_format = optarg;
			break;

		case 'h':
			usage(EXIT_SUCCESS);
			break;

		default:
			usage(EXIT_FAILURE);
		}
	}

	if (optind == argc) {
		usage(EXIT_FAILURE);
	}

	(void)mbdyn_make_salt(salt, sizeof(salt), salt_format);
	salt[STRLENOF(salt)] = '\0';

	if (asserted_cred) {
		if (strncmp(asserted_cred, "{CRYPT}", STRLENOF("{CRYPT}")) == 0) {
			asserted_cred += STRLENOF("{CRYPT}");
		} else {
			memcpy(asserted_cred_buf, asserted_cred, sizeof(asserted_cred_buf));
			asserted_cred_buf[STRLENOF(asserted_cred_buf)] = '\0';

			asserted_cred = crypt(asserted_cred_buf, salt);

			memcpy(asserted_cred_buf, asserted_cred, sizeof(asserted_cred_buf));
			asserted_cred[STRLENOF(asserted_cred_buf)] = '\0';
			asserted_cred = asserted_cred_buf;
		}
		asserted_cred_len = strlen(asserted_cred);
	}

	argv = &argv[optind];
	argc -= optind;

	for (; argc; argc--) {
		char	cred[35];
		char	*c;

		memcpy(cred, argv[0], sizeof(cred));
		cred[STRLENOF(cred)] = '\0';

		if (asserted_cred) {
			c = crypt(cred, asserted_cred);

			if (strcmp(c, asserted_cred) == 0) {
				printf("%s OK\n", c);
			} else {
				printf("%s ERR\n", c);
			}
		} else {
			c = crypt(cred, salt);
			printf("%s\n", c);
		}
	}

	return EXIT_SUCCESS;
#endif /* HAVE_CRYPT */	
}
