/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include "mbconfig.h"

#include <cstdlib>

#include <iostream>
#include <fstream>
#include "ac/getopt.h"

#include "myassert.h"
#include "mynewmem.h"

#include "aerodata.h"
#include "aerodc81.h"
#include "c81data.h"

static void
usage(int rc)
{
	std::cerr <<
"usage: c81test [<options>] <file>\n"
"\n"
"options:\n"
"\t-a <alpha>\t"	"dump coefficients for angle of attack <alpha>\n"
"\t-c\t\t"		"interpret <file> as traditional c81 format (default)\n"
"\t-C <coef>\t"		"only dump coefficient <coef> (must be cl, cd or cm)\n"
"\t-d {<dump>|-}\t"	"dump contents (to optional file <dump>, if given)\n"
"\t-f\t\t"		"interpret <file> as free format\n"
"\t-F\t\t"		"\"flip\" coefficients\n"
"\t-m <mach>\t"		"dump coefficients for Mach number <mach>\n"
"\t-n\t\t"		"interpret <file> as NREL format\n"
"\t-o\t\t"		"interpret <file> as fc511 format\n"
"\t-t <tolerance>\t"	"use <tolerance> to determine end of linearity\n"
		<< std::endl;

	exit(rc);
}

/*
 * uso: testc81 file alpha mach
 */
int
main(int argc, char *argv[])
{
	char		mode = 'x';
	char		*dump_fname = NULL;
	bool		dump(false);
	bool		flip(false);
	doublereal	alpha;
	doublereal	mach;
	enum {
		GOT_NONE	= 0x0000U,
		GOT_ALPHA	= 0x0001U,
		GOT_MACH	= 0x0002U,
		GOT_ALPHAMACH	= (GOT_ALPHA|GOT_MACH)
	};
	unsigned	got = GOT_NONE;
	char		coef = '\0';
	doublereal	tol = 1e-6;
	
	while (true) {
		int opt = getopt(argc, argv, "a:cC:d:fFm:not:");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'a':
			alpha = atof(optarg);
			got |= GOT_ALPHA;
			break;

		case 'c':
			mode = 'c';
			break;

		case 'C':
			if (strcasecmp(optarg, "cl") == 0) {
				coef = 'l';

			} else if (strcasecmp(optarg, "cd") == 0) {
				coef = 'd';

			} else if (strcasecmp(optarg, "cm") == 0) {
				coef = 'm';

			} else {
				std::cerr << "unknown coefficient "
					"\"" << optarg << "\"" << std::endl;
				usage(EXIT_FAILURE);
			}
			break;

		case 'd':
			dump = true;
			if (optarg && strcmp(optarg, "-") != 0) {
				dump_fname = optarg;
			}
			break;

		case 'f':
			mode = 'f';
			break;

		case 'F':
			flip = true;
			break;

		case 'm':
			mach = atof(optarg);
			got |= GOT_MACH;
			break;

		case 'n':
		case 'o':
			mode = char(opt);
			break;

		case 't':
			tol = atof(optarg);
			break;

		default:
			usage(EXIT_FAILURE);
		}
	}

	argc -= optind;
	argv += optind;

	/* C81 file name is mandatory */
	if (argc != 1) {
		usage(EXIT_FAILURE);
	}

	/* consistency checks... */
	if (dump_fname || dump) {
		if (got) {
			std::cerr << "\"dump\" incompatible "
				"with alpha/mach selection"
				<< std::endl;
			usage(EXIT_FAILURE);
		}

	} else {
		switch (got) {
		default:
		case GOT_NONE:
		case GOT_ALPHA:
		case GOT_MACH:
			if (coef == '\0') {
				std::cerr << "need to select a coefficient "
					"when alpha and mach "
					"are not selected altogether"
					<< std::endl;
				usage(EXIT_FAILURE);
			}
			break;

		case GOT_ALPHAMACH:
			break;
		}
	}

	/* read C81 database */
	std::ifstream in(argv[0]);
	if (!in) {
		std::cerr << "unable to open file '" << argv[0] << "'"
			<< std::endl;
		exit(EXIT_FAILURE);
	}

	c81_data databuf = {{'\0'}};
	c81_data* data = &databuf;

	int rc;

	int ff = 0;
	switch (mode) {
	case 'x':
		rc = c81_data_read(in, data, tol, &ff);
		break;

	case 'c':
		rc = c81_data_read(in, data, tol, 0);
		break;

	case 'f':
		rc = c81_data_read_free_format(in, data, tol);
		break;

	case 'n':
		rc = c81_data_nrel_read(in, data, tol);
		break;

	case 'o':
		rc = c81_data_fc511_read(in, data, tol);
		break;

	default:
		exit(EXIT_FAILURE);
	}

	in.close();

	if (rc) {
		std::cerr << "unable to read c81 data from file "
			"\"" << argv[0] << "\"" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (flip) {
		c81_data_flip(data);
	}

	if (dump_fname) {
		std::ofstream out(dump_fname);
		c81_data_write_free_format(out, data);
		out.close();
		c81_data_destroy(data);
		exit(EXIT_SUCCESS);

	} else if (dump) {
		c81_data_write_free_format(std::cout, data);
		c81_data_destroy(data);
		exit(EXIT_SUCCESS);
	}

	if ((got & GOT_ALPHAMACH) == GOT_ALPHAMACH) {
		switch (coef) {
		case 'l':
			std::cout << c81_data_get_coef(data->NML, data->ml,
				data->NAL, data->al, alpha, mach)
				<< std::endl;
			break;

		case 'd':
			std::cout << c81_data_get_coef(data->NMD, data->md,
				data->NAD, data->ad, alpha, mach)
				<< std::endl;
			break;

		case 'm':
			std::cout << c81_data_get_coef(data->NMM, data->mm,
				data->NAM, data->am, alpha, mach)
				<< std::endl;
			break;

		case '\0':
			std::cout 
				<< "alpha: " << alpha << std::endl
				<< "mach: " << mach << std::endl
				<< "cl: " << c81_data_get_coef(data->NML, data->ml,
					data->NAL, data->al, alpha, mach)
					<< std::endl
				<< "cd: " << c81_data_get_coef(data->NMD, data->md,
					data->NAD, data->ad, alpha, mach)
					<< std::endl
				<< "cm: " << c81_data_get_coef(data->NMM, data->mm,
					data->NAM, data->am, alpha, mach)
					<< std::endl;
			break;
		}

	} else if (got & GOT_ALPHAMACH) {
		int NM = 0, NA = 0;
		double *m = 0, *a = 0;

		switch (coef) {
		case 'l':
			NM = data->NML;
			NA = data->NAL;
			m = data->ml;
			a = data->al;
			break;

		case 'd':
			NM = data->NMD;
			NA = data->NAD;
			m = data->md;
			a = data->ad;
			break;

		case 'm':
			NM = data->NMM;
			NA = data->NAM;
			m = data->mm;
			a = data->am;
			break;

		case '\0':
			/* TODO */
			break;
		}

		if (m) {
			switch (got & GOT_ALPHAMACH) {
			case GOT_ALPHA:
				for (int i = 0; i < NM; i++) {
					std::cout
						<< m[i] << " "
						<< c81_data_get_coef(NM, m,
						NA, a, alpha, m[i])
						<< std::endl;
				}
				break;

			case GOT_MACH:
				for (int i = 0; i < NA; i++) {
					std::cout
						<< a[i] << " "
						<< c81_data_get_coef(NM, m,
						NA, a, a[i], mach)
						<< std::endl;
				}
				break;
			}
		}

	} else {
		std::cerr << "don't know what to do..." << std::endl;
		c81_data_destroy(data);
		usage(EXIT_FAILURE);
	}

	c81_data_destroy(data);

	return 0;
}
