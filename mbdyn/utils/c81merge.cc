/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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
#include <unistd.h>

#include <iostream>
#include <fstream>
#include "ac/getopt.h"

#include "myassert.h"
#include "mynewmem.h"

#include "aerodata.h"
#include "aerodc81.h"
#include "c81data.h"

static void
merge_vecs(int n, doublereal *v1, doublereal w1, doublereal *v2, doublereal w2, doublereal *v)
{
	for (int i = 0; i < n; i++) {
		v[i] = w1*v1[i] + w2*v2[i];
	}
}

static void
merge_data_from(int NA, int NM, int NM_to,
		doublereal wfrom, doublereal wto, doublereal *workv,
		doublereal *m_from, doublereal *m_to, doublereal *m_dst,
		doublereal *a_from, doublereal *a_to, doublereal *a_dst)
{
	int	mt = 0;

	for (int m = 0; m < NM; m++) {
		m_dst[m] = m_from[m];
	}

	for (int a = 0; a < NA; a++) {
		a_dst[a] = a_from[a];
	}

	for (int m = 0; m < NM; m++) {
		for (; mt < NM_to; mt++) {
			if (m_to[mt] >= m_dst[m]) {
				break;
			}
		}

		/* last, step back and maintain */
		if (mt == NM_to) {
			mt = NM_to - 1;
			merge_vecs(NA, &a_from[NA*(1 + m)], wfrom,
					&a_to[NA*(1 + mt)], wto,
					&a_dst[NA*(1 + m)]);
					
		/* first, or exact match */
		} else if (mt == 0 || m_to[mt] == m_dst[m]) {
			merge_vecs(NA, &a_from[NA*(1 + m)], wfrom,
					&a_to[NA*(1 + mt)], wto,
					&a_dst[NA*(1 + m)]);

		/* need interpolation */
		} else {
			doublereal dd = m_to[mt] - m_to[mt - 1];
			doublereal d1 = (m_to[mt] - m_dst[m])/dd;
			doublereal d2 = (m_dst[m] - m_to[mt - 1])/dd;
			merge_vecs(NA, &a_to[NA*(mt)], d1,
					&a_to[NA*(mt + 1)], d2,
					workv);
			merge_vecs(NA, &a_from[NA*(1 + m)], wfrom,
					workv, wto,
					&a_dst[NA*(1 + m)]);
		}
	}
}

static int
c81_data_merge_from(c81_data *data_from, doublereal wfrom, c81_data *data_to, doublereal wto, c81_data *data)
{
	*data = *data_from;

	data->ml = new doublereal[data->NML];
	data->al = new doublereal[data->NAL*(data->NML + 1)];

	data->md = new doublereal[data->NMD];
	data->ad = new doublereal[data->NAD*(data->NMD + 1)];

	data->mm = new doublereal[data->NMM];
	data->am = new doublereal[data->NAM*(data->NMM + 2)];

	int		NA, NM, NM_to;
	doublereal	*workv = &data->am[data->NAM*(data->NMM + 1)],
			*a_from, *a_to, *a_dst, *m_from, *m_to, *m_dst;

	a_from = data_from->al;
	a_to = data_to->al;
	a_dst = data->al;
	m_from = data_from->ml;
	m_to = data_to->ml;
	m_dst = data->ml;
	NA = data->NAL;
	NM = data->NML;
	NM_to = data_to->NML;
	merge_data_from(NA, NM, NM_to, wfrom, wto, workv, m_from, m_to, m_dst, a_from, a_to, a_dst);

	a_from = data_from->ad;
	a_to = data_to->ad;
	a_dst = data->ad;
	m_from = data_from->md;
	m_to = data_to->md;
	m_dst = data->md;
	NA = data->NAD;
	NM = data->NMD;
	NM_to = data_to->NMD;
	merge_data_from(NA, NM, NM_to, wfrom, wto, workv, m_from, m_to, m_dst, a_from, a_to, a_dst);

	a_from = data_from->am;
	a_to = data_to->am;
	a_dst = data->am;
	m_from = data_from->mm;
	m_to = data_to->mm;
	m_dst = data->mm;
	NA = data->NAM;
	NM = data->NMM;
	NM_to = data_to->NMM;
	merge_data_from(NA, NM, NM_to, wfrom, wto, workv, m_from, m_to, m_dst, a_from, a_to, a_dst);

	// c81_do_data_stall(data);

	return 0;
}

static void
merge_data_to(int NA, int NM, int NM_from,
		doublereal wfrom, doublereal wto, doublereal *workv,
		doublereal *m_from, doublereal *m_to, doublereal *m_dst,
		doublereal *a_from, doublereal *a_to, doublereal *a_dst)
{
	int	mf = 0;

	for (int m = 0; m < NM; m++) {
		m_dst[m] = m_to[m];
	}

	for (int a = 0; a < NA; a++) {
		a_dst[a] = a_to[a];
	}

	for (int m = 0; m < NM; m++) {
		for (; mf < NM_from; mf++) {
			if (m_from[mf] >= m_dst[m]) {
				break;
			}
		}

		/* last, step back and maintain */
		if (mf == NM_from) {
			mf = NM_from - 1;
			merge_vecs(NA, &a_from[NA*(1 + mf)], wfrom,
					&a_to[NA*(1 + m)], wto,
					&a_dst[NA*(1 + m)]);

		/* first, or exact match */
		} else if (mf == 0 || m_from[mf] == m_dst[m]) {
			merge_vecs(NA, &a_from[NA*(1 + mf)], wfrom,
					&a_to[NA*(1 + m)], wto,
					&a_dst[NA*(1 + m)]);

		/* need interpolation */
		} else {
			doublereal dd = m_from[mf] - m_from[mf - 1];
			doublereal d1 = (m_from[mf] - m_dst[m])/dd;
			doublereal d2 = (m_dst[m] - m_from[mf - 1])/dd;
			merge_vecs(NA, &a_from[NA*(mf)], d1,
					&a_from[NA*(mf + 1)], d2,
					workv);
			merge_vecs(NA, workv, wfrom,
					&a_to[NA*(1 + m)], wto,
					&a_dst[NA*(1 + m)]);
		}
	}
}

static int
c81_data_merge_to(c81_data *data_from, doublereal wfrom, c81_data *data_to, doublereal wto, c81_data *data)
{
	*data = *data_to;

	data->ml = new doublereal[data->NML];
	data->al = new doublereal[data->NAL*(data->NML + 1)];

	data->md = new doublereal[data->NMD];
	data->ad = new doublereal[data->NAD*(data->NMD + 1)];

	data->mm = new doublereal[data->NMM];
	data->am = new doublereal[data->NAM*(data->NMM + 2)];

	int		NA, NM, NM_from;
	doublereal	*workv = &data->am[data->NAM*(data->NMM + 1)],
			*a_from, *a_to, *a_dst, *m_from, *m_to, *m_dst;

	a_from = data_from->al;
	a_to = data_to->al;
	a_dst = data->al;
	m_from = data_from->ml;
	m_to = data_to->ml;
	m_dst = data->ml;
	NA = data->NAL;
	NM = data->NML;
	NM_from = data_from->NML;
	merge_data_to(NA, NM, NM_from, wfrom, wto, workv, m_from, m_to, m_dst, a_from, a_to, a_dst);

	a_from = data_from->ad;
	a_to = data_to->ad;
	a_dst = data->ad;
	m_from = data_from->md;
	m_to = data_to->md;
	m_dst = data->md;
	NA = data->NAD;
	NM = data->NMD;
	NM_from = data_from->NMD;
	merge_data_to(NA, NM, NM_from, wfrom, wto, workv, m_from, m_to, m_dst, a_from, a_to, a_dst);

	a_from = data_from->am;
	a_to = data_to->am;
	a_dst = data->am;
	m_from = data_from->mm;
	m_to = data_to->mm;
	m_dst = data->mm;
	NA = data->NAM;
	NM = data->NMM;
	NM_from = data_from->NMM;
	merge_data_to(NA, NM, NM_from, wfrom, wto, workv, m_from, m_to, m_dst, a_from, a_to, a_dst);

	// c81_do_data_stall(data);

	return 0;
}

static int
c81_data_merge_both(c81_data *data_from, doublereal wfrom, c81_data *data_to, doublereal wto, c81_data *data)
{
	return -1;
}

/*
 * uso: 
 */
int
main(int argc, char *argv[])
{
	c81_data	data_from,
			data_to,
			data;
	char		*name_from = 0,
			*name_to = 0,
			*name = 0,
			header[31] = { '\0' },
			*next;
	doublereal	dFrom = -1.,
			dTo;
	int		rc = EXIT_SUCCESS;
	doublereal	tol = 1e-6;
	
	enum {
		MODE_UNDEFINED, MODE_FROM, MODE_TO, MODE_BOTH
	} mode = MODE_UNDEFINED;

	while (1) {
		int	opt = getopt(argc, argv, "f:hH:m:o:s:t:");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'f':
			if (name_from != 0) {
				silent_cerr("-f already set to \"" << name_from << "\"" << std::endl);
			}
			name_from = optarg;
			break;

		case 'H':
			if (header[0] != '\0') {
				silent_cerr("-h already set to \"" << header << "\"" << std::endl);
			}
			memcpy(header, optarg, sizeof(header));
			header[STRLENOF(header)] = '\0';
			break;

		case 'm':
			switch (mode) {
			case MODE_FROM:
				silent_cerr("-m already set to \"FROM\"" << std::endl);
				break;

			case MODE_TO:
				silent_cerr("-m already set to \"TO\"" << std::endl);
				break;

			case MODE_BOTH:
				silent_cerr("-m already set to \"BOTH\"" << std::endl);
				break;

			case MODE_UNDEFINED:
				break;
			}

			if (strcasecmp(optarg, "from") == 0) {
				mode = MODE_FROM;
				
			} else if (strcasecmp(optarg, "to") == 0) {
				mode = MODE_TO;
				
			} else if (strcasecmp(optarg, "both") == 0) {
				mode = MODE_BOTH;
				
			} else {
				silent_cerr("unknown mode \"" << optarg << "\"" <<std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			break;

		case 'o':
			if (name != 0) {
				silent_cerr("-o already set to \"" << name << "\"" << std::endl);
			}
			name = optarg;
			break;

		case 's':
			if (dFrom != -1.) {
				silent_cerr("-s already set to " << dFrom << std::endl);
			}
			dFrom = strtod(optarg, &next);
			if (next == optarg || next[0] != '\0') {
				silent_cerr("illegal value \"" << optarg << "\" for -s option" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			if (dFrom < 0. || dFrom > 1.) {
				silent_cerr("-s " << optarg << " is out of bounds (0.,1.)" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			break;

		case 't':
			if (name_to != 0) {
				silent_cerr("-t already set to \"" << name_to << "\"" << std::endl);
			}
			name_to = optarg;
			break;

		default:
			rc = EXIT_FAILURE;

		case 'h':
			silent_cout(
"\n"
"c81merge: merges two c81 files in a given proportion\n"
"\n"
"\t-f first.c81\t"	"first c81 data set\n"
"\t-h\t\t"		"this message\n"
"\t-H header\t"		"the header (max 30 chars)\n"
"\t-m {from|to|both}\t"	"use alpha/mach of \"from\", \"to\" or both\n"
"\t-o out.c81\t"	"output file name (stdout if missing)\n"
"\t-s s in [0,1]\t"	"real number; result = s * first + (1 - s) * second\n"
"\t-t second.c81\t"	"second c81 data set\n" 
"\n"
					);
			exit(rc);
		}
	}

	if (name_from == 0) {
		silent_cerr("missing required -f \"from\" parameter" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (name_to == 0) {
		silent_cerr("missing required -t \"to\" parameter" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (dFrom == -1.) {
		silent_cerr("missing required -s \"fraction\" parameter" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (mode == MODE_UNDEFINED) {
		mode = MODE_BOTH;
	}

	dTo = 1. - dFrom;

	std::ifstream in;
	
	/* from */
	in.open(name_from);
	if (!in) {
		silent_cerr("unable to open file \"" << name_from << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	int ff_from = 0;
	if (c81_data_read(in, &data_from, tol, &ff_from)) {
		silent_cerr("unable to read c81 data from file \"" << name_from << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	in.close();

	/* to */
	in.open(name_to);
	if (!in) {
		silent_cerr("unable to open file \"" << name_to << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	int ff_to = 0;
	if (c81_data_read(in, &data_to, tol, &ff_to)) {
		silent_cerr("unable to read c81 data from file \"" << name_to << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	in.close();

	/* consistency: at the moment, the two data sets must have exactly the same AoA pattern */
	if (data_from.NAL != data_to.NAL) {
		silent_cerr("number of AoA values for Cl differ ("
				<< data_from.NAL << " vs. " << data_to.NAL
				<< ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (int i = 0; i < data_from.NAL; i++) {
		if (data_from.al[i] != data_to.al[i]) {
			silent_cerr("AoA value " << i << " for Cl differs ("
					<< data_from.al[i] << " vs. " 
					<< data_to.al[i] << ")"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (data_from.NAD != data_to.NAD) {
		silent_cerr("number of AoA values for Cd differ ("
				<< data_from.NAD << " vs. " << data_to.NAD
				<< ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (int i = 0; i < data_from.NAL; i++) {
		if (data_from.ad[i] != data_to.ad[i]) {
			silent_cerr("AoA value " << i << " for Cd differs ("
					<< data_from.ad[i] << " vs. " 
					<< data_to.ad[i] << ")"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if (data_from.NAM != data_to.NAM) {
		silent_cerr("number of AoA values for Cm differ ("
				<< data_from.NAM << " vs. " << data_to.NAM
				<< ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (int i = 0; i < data_from.NAL; i++) {
		if (data_from.am[i] != data_to.am[i]) {
			silent_cerr("AoA value " << i << " for Cm differs ("
					<< data_from.am[i] << " vs. " 
					<< data_to.am[i] << ")"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	/* merge */
	switch (mode) {
	case MODE_FROM:
		rc = c81_data_merge_from(&data_from, dFrom, &data_to, dTo, &data);
		break;

	case MODE_TO:
		rc = c81_data_merge_to(&data_from, dFrom, &data_to, dTo, &data);
		break;

	case MODE_BOTH:
		rc = c81_data_merge_both(&data_from, dFrom, &data_to, dTo, &data);
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	memcpy(data.header, header, sizeof(header));

	if (rc) {
		silent_cerr("data merge failed" << std::endl);
		return rc;
	}

	if (name) {
		std::ofstream fout(name);
		if (!fout) {
			silent_cerr("unable to open file \"" << name << "\"" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		rc = c81_data_write(fout, &data);

	} else {
		rc = c81_data_write(std::cout, &data);
	}

	if (rc) {
		silent_cerr("output failed" << std::endl);
	}
	
	return rc;
}
