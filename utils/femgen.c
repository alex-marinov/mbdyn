/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "mbconfig.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ac/f2c.h"

extern int __FC_DECL__(femgen)(char outname[73], int32_t *iimd, int32_t *iimv, int32_t *idxm);

static void
usage(FILE *outf, int rc)
{
	fprintf(outf,
"usage: femgen [-hd] [-m {cxyz}] [[-o] <outfile>]\n"
"\t-d\t\tno initial modal displacements/velocities\n"
"\t-h\t\tthis message\n"
"\t-m <idx>\tindex of \"mass\" component to be used\n"
"\t\t\t('x', 'y', 'z'; 'c' to check consistency)\n"
"\t-o <outfile>\toutput file name (up to 72 characters long)\n"
		);
	exit(rc);

}

int
main(int argc, char *argv[])
{
	char outname[73] = { ' ' };

	int32_t iimd = 1, iimv = 1, idxm = 0;

	for (;;) {
		int opt = getopt(argc, argv, "dhm:o:");
		if (opt == -1) {
			break;
		}

		switch (opt) {
		case 'd':
			iimd = 0;
			iimv = 0;
			break;

		case '?':
		case 'h':
			usage(stdout, EXIT_SUCCESS);

		case 'm':
			if (optarg[1] != '\0') {
				usage(stderr, EXIT_FAILURE);
			}

			switch (optarg[0]) {
			case 'c':
				idxm = -1;
				break;

			case 'x':
				idxm = 1;
				break;

			case 'y':
				idxm = 2;
				break;

			case 'z':
				idxm = 3;
				break;

			default:
				fprintf(stderr, "femgen: unhandled parameter '%c' for option '-m'\n", optarg[0]);
				usage(stderr, EXIT_FAILURE);
			}
			break;

		case 'o': {
			size_t len = strlen(optarg);
			if (len >= sizeof(outname)) {
				fprintf(stderr, "femgen: output file name '%s' too long; trim to 72 bytes or less\n", optarg);
				exit(EXIT_FAILURE);
			}

			strcpy(outname, optarg);
			} break;

		default:
			fprintf(stderr, "femgen: unhandled option '%c'\n", opt);
			usage(stderr, EXIT_FAILURE);
		}
	}

	if (optind < argc) {
		if (outname[0] != ' ') {
			fprintf(stderr, "femgen: output file name already set using '-o' option\n");
			exit(EXIT_FAILURE);

		} else {
			size_t len = strlen(argv[optind]);
			if (len >= sizeof(outname)) {
				fprintf(stderr, "femgen: output file name '%s' too long; trim to 72 bytes or less\n", argv[optind]);
				exit(EXIT_FAILURE);
			}

			strcpy(outname, argv[optind]);
			optind++;
		}
	}

	if (optind < argc) {
		fprintf(stderr, "femgen: extra args ignored\n");
	}

	return __FC_DECL__(femgen)(outname, &iimd, &iimv, &idxm);
}

