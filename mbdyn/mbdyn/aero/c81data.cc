/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "ac/f2c.h"
#include "myassert.h"

extern "C" {
#include "aerodc81.h"
}

#include "c81data.h"

/*
 * header NML,NAL,NMD,NAD,NMM,NAM	A30,6I2
 * ML(1),....,ML(NML)			7x,9F7.0   eventualmente su piu' righe
 * AL(1)  CL(1,1),....,CL(1,NML)	10F7.0/(7x,9F7.0)
 * :         :     :       :
 * AL(NAL)CL(NAL,1),....,CL(NAL,NML)	10F7.0/(7x,9F7.0)
 * AD(1)  CD(1,1),....,CD(1,NMD)	10F7.0/(7x,9F7.0)
 * :         :     :       :
 * AD(NAD)CD(NAD,1),....,CD(NAD,NMD)	10F7.0/(7x,9F7.0)
 * AM(1)  CM(1,1),....,CL(1,NMM)	10F7.0/(7x,9F7.0)
 * :         :     :       :
 * AM(NAM)CM(NAM,1),....,CL(NAM,NMM)	10F7.0/(7x,9F7.0)
 */

static int
do_stall(int NM, int NA, doublereal *a, doublereal *stall, const doublereal dcltol);

static int
get_int(const char *const from, int &i)
{
#ifdef HAVE_STRTOL
	char *endptr = NULL;
	errno = 0;
	i = strtol(from, &endptr, 10);
	int save_errno = errno;
	if (endptr != NULL && endptr[0] != '\0' && !isspace(endptr[0])) {
		return -1;

	} else if (save_errno == ERANGE) {
		silent_cerr("c81data: warning, int "
			<< std::string(from, endptr - from)
			<< " overflows" << std::endl);
		return -1;
	}
#else /* !HAVE_STRTOL */
   	i = atoi(buf);
#endif /* !HAVE_STRTOL */
	return 0;
}

#if 0	// unused so far
static int
get_long(const char *const from, long &l)
{
#ifdef HAVE_STRTOL
	char *endptr = NULL;
	errno = 0;
	l = strtol(from, &endptr, 10);
	int save_errno = errno;
	if (endptr != NULL && endptr[0] != '\0' && !isspace(endptr[0])) {
		return -1;

	} else if (save_errno == ERANGE) {
		silent_cerr("c81data: warning, int "
			<< std::string(from, endptr - from)
			<< " overflows" << std::endl);
		return -1;
	}
#else /* !HAVE_STRTOL */
   	l = atol(buf);
#endif /* !HAVE_STRTOL */
	return 0;
}
#endif

static int
get_double(const char *const from, doublereal &d)
{
#ifdef HAVE_STRTOD
	char *endptr = NULL;
	errno = 0;
	d = strtod(from, &endptr);
	int save_errno = errno;
	if (endptr != NULL && endptr[0] != '\0' && !isspace(endptr[0])) {
		return -1;

	} else if (save_errno == ERANGE) {
		silent_cerr("c81data: warning, double "
			<< std::string(from, endptr - from)
			<< " overflows" << std::endl);
		return -1;
	}
#else /* !HAVE_STRTOD */
   	d = atof(buf);
#endif /* !HAVE_STRTOD */
	return 0;
}

static int
get_c81_vec(std::istream& in, doublereal* v, int ncols)
{
	char	buf[128];

   	if (!in || v == NULL || ncols < 1) {
      		return -1;
   	}

   	for (int c = 0; c < ncols; c++) {
		char	tmp[8];
		int	i = c%9;

		if (i == 0) {
			in.getline(buf, sizeof(buf));
		}

		/* always skip first */
		memcpy(tmp, &buf[7*(i + 1)], 7);
		tmp[7] = '\0';

		if (get_double(tmp, v[c])) {
			return -1;
		}
   	}

   	return 0;
}

static int
check_vec(doublereal* v, int n)
{
	for (int i = 1; i < n; i++) {
		if (v[i] <= v[i - 1]) {
			return i;
		}
	}

	return 0;
}

static int
get_c81_mat(std::istream& in, doublereal* m, int nrows, int ncols)
{
	char	buf[128];

   	if (!in || m == NULL || nrows < 1 || ncols < 1) {
      		return -1;
   	}

   	for (int r = 0; r < nrows; r++) {
		int	row = 0;
		int	offset = 0;

      		for (int c = 0; c < ncols; c++) {
			char	tmp[8];
			int	i = (c - offset)%(9 - offset + 1);

			if (i == 0) {
				if (++row == 2) {
					offset = 1;
				}
				in.getline(buf, sizeof(buf));
			}

			/* skip first except first time */
			memcpy(tmp, &buf[7*(i + offset)], 7);
			tmp[7] = '\0';

			if (get_double(tmp, m[r + nrows*c])) {
				return -1;
			}
      		}
   	}

   	return 0;
}

static int
put_c81_row(std::ostream& out, doublereal* v, int dim, int ncols, int first = 0)
{
   	int start = 0;
   	const int N = 9;

   	if (first) {
      		out << std::setw(7) << v[0];
      		start = dim;
      		ncols--;
   	} else {
      		out << std::setw(7) << "";
   	}

   	for (int i = 0; i < (ncols-1)/N; i++) {
      		for (int j = 0; j < N; j++) {
	 		out << std::setw(7) << v[start+dim*j];
     	 	}
      		out << std::endl << std::setw(7) << "";
      		start += dim*N;
   	}

   	for (int j = 0; j < (ncols-1)%N+1; j++) {
      		out << std::setw(7) << v[start+dim*j];
   	}
   	out << std::endl;

   	return 0;
}

static int
put_c81_vec(std::ostream& out, doublereal* v, int nrows)
{
   	if (!out || v == NULL || nrows < 1) {
      		return -1;
   	}

   	put_c81_row(out, v, 1, nrows);

   	return 0;
}

static int
put_c81_mat(std::ostream& out, doublereal* m, int nrows, int ncols)
{
   	if (!out || m == NULL || nrows < 1 || ncols < 1) {
      		return -1;
   	}

   	for (int i = 0; i < nrows; i++) {
      		put_c81_row(out, m+i, nrows, ncols, 1);
   	}

   	return 0;
}

extern "C" void
c81_data_destroy(c81_data* data)
{
	delete[] data->ml;
	delete[] data->al;

	delete[] data->md;
	delete[] data->ad;

	delete[] data->mm;
	delete[] data->am;

	delete[] data->stall;
	delete[] data->mstall;
}

extern "C" int
c81_data_read(std::istream& in, c81_data* data, const doublereal dcltol, int *ff)
{
   	char buf[BUFSIZ];	// 81 should suffice

	if (ff) {
		*ff = 0;
	}

   	/* header */
   	in.getline(buf, sizeof(buf));
	if (strncasecmp(buf, "# FREE FORMAT", STRLENOF("# FREE FORMAT")) == 0) {
		if (ff) {
			*ff = 1;
		}

		return c81_data_read_free_format(in, data, dcltol);
	}

   	buf[42] = '\0';
   	if (get_int(&buf[40], data->NAM)) {
		return -1;
	}

   	buf[40] = '\0';
   	if (get_int(&buf[38], data->NMM)) {
		return -1;
	}

   	buf[38] = '\0';
   	if (get_int(&buf[36], data->NAD)) {
		return -1;
	}

   	buf[36] = '\0';
   	if (get_int(&buf[34], data->NMD)) {
		return -1;
	}

   	buf[34] = '\0';
   	if (get_int(&buf[32], data->NAL)) {
		return -1;
	}

   	buf[32] = '\0';
   	if (get_int(&buf[30], data->NML)) {
		return -1;
	}

	if (data->NAM <= 0 || data->NMM <= 0
			|| data->NAD <= 0 || data->NMD <= 0
			|| data->NAL <= 0 || data->NML <= 0) {
		return -1;
	}

   	buf[30] = '\0';
   	strncpy(data->header, buf, 30);
   	data->header[30] = '\0';

   	/* lift */
   	data->ml = new doublereal[data->NML];
   	get_c81_vec(in, data->ml, data->NML);
	if (check_vec(data->ml, data->NML)) {
		return -1;
	}

   	data->al = new doublereal[(data->NML+1)*data->NAL];
   	get_c81_mat(in, data->al, data->NAL, data->NML+1);
	if (check_vec(data->al, data->NAL)) {
		return -1;
	}

   	/* drag */
   	data->md = new doublereal[data->NMD];
   	get_c81_vec(in, data->md, data->NMD);
	if (check_vec(data->md, data->NMD)) {
		return -1;
	}

   	data->ad = new doublereal[(data->NMD+1)*data->NAD];
   	get_c81_mat(in, data->ad, data->NAD, data->NMD+1);
	if (check_vec(data->ad, data->NAD)) {
		return -1;
	}

   	/* moment */
   	data->mm = new doublereal[data->NMM];
   	get_c81_vec(in, data->mm, data->NMM);
	if (check_vec(data->mm, data->NMM)) {
		return -1;
	}

   	data->am = new doublereal[(data->NMM+1)*data->NAM];
   	get_c81_mat(in, data->am, data->NAM, data->NMM+1);
	if (check_vec(data->am, data->NAM)) {
		return -1;
	}

	/* FIXME: maybe this is not the best place */
	c81_data_do_stall(data, dcltol);

   	return 0;
}

extern "C" int
c81_data_merge(
	unsigned ndata,
	const c81_data **data,
	const doublereal *upper_bounds,
	doublereal dCsi,
	doublereal dcltol,
	c81_data *i_data)
{
	ASSERT(ndata > 0);
	ASSERT(data != 0);
	ASSERT(upper_bounds != 0);
	ASSERT(dCsi >= -1.);
	ASSERT(dCsi <= 1.);
	ASSERT(dcltol > 0.);
	ASSERT(i_data != 0);

	if (dCsi < upper_bounds[0]) {
		silent_cerr("cannot find C81 data lower bound for point xi=" << dCsi << std::endl);
		return -1;
	}

	int to = 0;
	for (unsigned i = 1; i < ndata; i++) {
		if (upper_bounds[i] > dCsi) {
			to = i;
			break;
		}
	}
	if (to == 0) {
		silent_cerr("cannot find C81 data upper bound for point xi=" << dCsi << std::endl);
		return -1;
	}

	int from = to - 1;
	if (unsigned(from) == ndata) {
		silent_cerr("cannot find C81 data upper bound for point xi=" << dCsi << std::endl);
		return -1;
	}

	/* we need to interpolate between data[from]
	 * and data[to] */

	/* we only accept homogeneous data sources,
	 * i.e. same Mach and alpha patterns */
	if (data[from]->NML != data[to]->NML) {
		silent_cerr("number of Mach points for Cl between airfoils "
			<< from << " (" << data[from]->NML << ") and "
			<< to << " (" << data[to]->NML << ") do not match"
			<< std::endl);
		return -1;
	}

	if (data[from]->NAL != data[to]->NAL) {
		silent_cerr("number of AoA points for Cl between airfoils "
			<< from << " (" << data[from]->NAL << ") and "
			<< to << " (" << data[to]->NAL << ") do not match"
			<< std::endl);
		return -1;
	}

	if (data[from]->NMD != data[to]->NMD) {
		silent_cerr("number of Mach points for Cd between airfoils "
			<< from << " (" << data[from]->NMD << ") and "
			<< to << " (" << data[to]->NMD << ") do not match"
			<< std::endl);
		return -1;
	}

	if (data[from]->NAD != data[to]->NAD) {
		silent_cerr("number of AoA points for Cd between airfoils "
			<< from << " (" << data[from]->NAD << ") and "
			<< to << " (" << data[to]->NAD << ") do not match"
			<< std::endl);
		return -1;
	}

	if (data[from]->NMM != data[to]->NMM) {
		silent_cerr("number of Mach points for Cm between airfoils "
			<< from << " (" << data[from]->NMM << ") and "
			<< to << " (" << data[to]->NMM << ") do not match"
			<< std::endl);
		return -1;
	}

	if (data[from]->NAM != data[to]->NAM) {
		silent_cerr("number of AoA points for Cm between airfoils "
			<< from << " (" << data[from]->NAM << ") and "
			<< to << " (" << data[to]->NAM << ") do not match"
			<< std::endl);
		return -1;
	}

	for (int i = 0; i < data[from]->NML; i++) {
		if (data[from]->ml[i] != data[to]->ml[i]) {
			silent_cerr("Mach point " << i << "for Cl of airfoils "
				<< from << " (" << data[from]->ml[i] << ") and "
				<< to << " (" << data[to]->ml[i] << ") differs"
				<< std::endl);
			return -1;
		}
	}

	for (int i = 0; i < data[from]->NAL; i++) {
		if (data[from]->al[i] != data[to]->al[i]) {
			silent_cerr("AoA point " << i << "for Cl of airfoils "
				<< from << " (" << data[from]->al[i] << ") and "
				<< to << " (" << data[to]->al[i] << ") differs"
				<< std::endl);
			return -1;
		}
	}

	for (int i = 0; i < data[from]->NMD; i++) {
		if (data[from]->md[i] != data[to]->md[i]) {
			silent_cerr("Mach point " << i << "for Cd of airfoils "
				<< from << " (" << data[from]->md[i] << ") and "
				<< to << " (" << data[to]->md[i] << ") differs"
				<< std::endl);
			return -1;
		}
	}

	for (int i = 0; i < data[from]->NAD; i++) {
		if (data[from]->ad[i] != data[to]->ad[i]) {
			silent_cerr("AoA point " << i << "for Cd of airfoils "
				<< from << " (" << data[from]->ad[i] << ") and "
				<< to << " (" << data[to]->ad[i] << ") differs"
				<< std::endl);
			return -1;
		}
	}

	for (int i = 0; i < data[from]->NMM; i++) {
		if (data[from]->mm[i] != data[to]->mm[i]) {
			silent_cerr("Mach point " << i << "for Cm of airfoils "
				<< from << " (" << data[from]->mm[i] << ") and "
				<< to << " (" << data[to]->mm[i] << ") differs"
				<< std::endl);
			return -1;
		}
	}

	for (int i = 0; i < data[from]->NAM; i++) {
		if (data[from]->am[i] != data[to]->am[i]) {
			silent_cerr("AoA point " << i << "for Cm of airfoils "
				<< from << " (" << data[from]->am[i] << ") and "
				<< to << " (" << data[to]->am[i] << ") differs"
				<< std::endl);
			return -1;
		}
	}

	snprintf(i_data->header, sizeof(i_data->header),
		"interpolated (\"%s\"[%e]->\"%s\"[%e]: [%e])",
		data[from]->header, upper_bounds[from],
		data[to]->header, upper_bounds[to], dCsi);

	i_data->NML = data[from]->NML;
	i_data->NAL = data[from]->NAL;
	i_data->NMD = data[from]->NMD;
	i_data->NAD = data[from]->NAD;
	i_data->NMM = data[from]->NMM;
	i_data->NAM = data[from]->NAM;

	doublereal dw = upper_bounds[to] - upper_bounds[from];
	doublereal dw_from = (upper_bounds[to] - dCsi)/dw;
	doublereal dw_to = (dCsi - upper_bounds[from])/dw;

	/* lift */
	i_data->ml = new doublereal[i_data->NML];
	for (int i = 0; i < i_data->NML; i++) {
		i_data->ml[i] = data[from]->ml[i];
	}
	i_data->al = new doublereal[(i_data->NML + 1)*i_data->NAL];
	for (int i = 0; i < i_data->NAL; i++) {
		i_data->al[i] = data[from]->al[i];
	}
	for (int i = i_data->NAL; i < (i_data->NML + 1)*i_data->NAL; i++) {
		i_data->al[i] = dw_from*data[from]->al[i] + dw_to*data[to]->al[i];
	}

	/* drag */
	i_data->md = new doublereal[i_data->NMD];
	for (int i = 0; i < i_data->NMD; i++) {
		i_data->md[i] = data[from]->md[i];
	}
	i_data->ad = new doublereal[(i_data->NMD + 1)*i_data->NAD];
	for (int i = 0; i < i_data->NAD; i++) {
		i_data->ad[i] = data[from]->ad[i];
	}
	for (int i = i_data->NAD; i < (i_data->NMD + 1)*i_data->NAD; i++) {
		i_data->ad[i] = dw_from*data[from]->ad[i] + dw_to*data[to]->ad[i];
	}

	/* moment */
	i_data->mm = new doublereal[i_data->NMM];
	for (int i = 0; i < i_data->NMM; i++) {
		i_data->mm[i] = data[from]->mm[i];
	}
	i_data->am = new doublereal[(i_data->NMM + 1)*i_data->NAM];
	for (int i = 0; i < i_data->NAM; i++) {
		i_data->am[i] = data[from]->am[i];
	}
	for (int i = i_data->NAM; i < (i_data->NMM + 1)*i_data->NAM; i++) {
		i_data->am[i] = dw_from*data[from]->am[i] + dw_to*data[to]->am[i];
	}

	// FIXME: maybe this is not the best place
	c81_data_do_stall(i_data, dcltol);

	return 0;
}

extern "C" int
c81_data_write(std::ostream& out, c81_data* data)
{
	if (data == 0) {
		return -1;
	}

	std::ios::fmtflags tmpflags;
	tmpflags = out.flags(std::ios::left);
   	out << std::setw(30) << data->header;
	out.flags(tmpflags);
	out
     		<< std::setw(2) << data->NML
     		<< std::setw(2) << data->NAL
     		<< std::setw(2) << data->NMD
     		<< std::setw(2) << data->NAD
     		<< std::setw(2) << data->NMM
     		<< std::setw(2) << data->NAM
		<< std::endl;

   	put_c81_vec(out, data->ml, data->NML);
   	put_c81_mat(out, data->al, data->NAL, data->NML+1);

   	put_c81_vec(out, data->md, data->NMD);
   	put_c81_mat(out, data->ad, data->NAD, data->NMD+1);

   	put_c81_vec(out, data->mm, data->NMM);
   	put_c81_mat(out, data->am, data->NAM, data->NMM+1);

   	return 0;
}

extern "C" int
read_fc511_row(std::istream& in, doublereal *d, int NC)
{
	int	r;

	for (r = 0; r < NC/8; r++) {
		char	buf[81];

		in.getline(buf, sizeof(buf));

		for (int c = 7; c >= 0; c--) {
			buf[10*(c + 1)] = '\0';
		   	if (get_double(&buf[10*c], d[8*r + c])) {
				return -1;
			}
		}
	}

	if (NC%8) {
		char	buf[81];

		in.getline(buf, sizeof(buf));

		for (int c = NC%8 - 1; c >= 0; c--) {
			buf[10*(c + 1)] = '\0';
		   	if (get_double(&buf[10*c], d[8*r + c])) {
				return -1;
			}
		}
	}

	return 0;
}

extern "C" int
read_fc511_mat(std::istream& in, doublereal *d, int NR, int NC)
{
	for (int i = 0; i < NR; i++) {
		int	r;

		for (r = 0; r < NC/8; r++) {
			char	buf[81];

			in.getline(buf, sizeof(buf));

			for (int c = 7; c >= 0; c--) {
				buf[10*(c + 1)] = '\0';
			   	if (get_double(&buf[10*c], d[NR*(8*r + c) + i])) {
					return -1;
				}
			}
		}

		if (NC%8) {
			char	buf[81];

			in.getline(buf, sizeof(buf));

			for (int c = NC%8 - 1; c >= 0; c--) {
				buf[10*(c + 1)] = '\0';
			   	if (get_double(&buf[10*c], d[NR*(8*r + c) + i])) {
					return -1;
				}
			}
		}
	}

	return 0;
}

extern "C" int
read_fc511_block(std::istream& in, int &NA, int &NM, doublereal *&da, doublereal *&dm)
{
   	char buf[128];	// 81 should suffice; let's make it 128

   	in.getline(buf, sizeof(buf));

   	buf[10] = '\0';
   	if (get_int(&buf[5], NM)) {
		return -1;
	}

   	buf[5] = '\0';
   	if (get_int(&buf[0], NA)) {
		return -1;
	}

	if (NA <= 0 || NM <= 0) {
		return -1;
	}

	dm = new doublereal[NM];
	doublereal *tmpda = new doublereal[NA*(NM + 1)];

	if (read_fc511_row(in, tmpda, NA)) {
		return -1;
	}

	if (read_fc511_row(in, dm, NM)) {
		return -1;
	}

	if (read_fc511_mat(in, &tmpda[NA], NA, NM)) {
		return -1;
	}

	int pos;
   	in.getline(buf, sizeof(buf));
   	if (get_int(buf, pos)) {
		return -1;
	}

	doublereal *posda = new doublereal[2*pos];

	if (read_fc511_row(in, posda, pos)) {
		return -1;
	}

	if (read_fc511_row(in, &posda[pos], pos)) {
		return -1;
	}

	int neg;
   	in.getline(buf, sizeof(buf));
   	if (get_int(buf, neg)) {
		return -1;
	}

	int realNA = (neg + NA + pos);
	da = new doublereal[realNA*(NM + 1)];

	if (read_fc511_row(in, da, neg)) {
		return -1;
	}

	if (read_fc511_row(in, &da[realNA], neg)) {
		return -1;
	}

	for (int c = 0; c < 2; c++) {
		for (int r = 0; r < pos; r++) {
			da[realNA*c + neg + NA + r] = posda[pos*c + r];
		}
	}

	for (int c = 2; c < NM + 1; c++) {
		for (int r = 0; r < neg; r++) {
			da[realNA*c + r] = da[realNA + r];
		}
		for (int r = 0; r < pos; r++) {
			da[realNA*c + neg + NA + r] = posda[pos + r];
		}
	}

	for (int c = 0; c < (NM + 1); c++) {
		for (int r = 0; r < NA; r++) {
			da[realNA*c + neg + r] = tmpda[NA*c + r];
		}
	}

	delete[] tmpda;
	delete[] posda;

	NA = realNA;

	return 0;
}

extern "C" int
c81_data_fc511_read(std::istream& in, c81_data* data, const doublereal dcltol)
{
   	char buf[128];	// 81 should suffice; let's make it 128

   	/* header */
   	in.getline(buf, sizeof(buf));

	memcpy(data->header, buf, sizeof(data->header));
	data->header[STRLENOF(data->header)] = '\0';

	char	*p;

   	in.getline(buf, sizeof(buf));
	p = buf;

	while (isspace(p[0])) {
		p++;
	}

	if (!p[0] || strncasecmp(p, "cl", 2) != 0 || (p[2] != '\0' && !isspace(p[2]))) {
		return -1;
	}

	if (read_fc511_block(in, data->NAL, data->NML, data->al, data->ml)) {
		return -1;
	}

   	in.getline(buf, sizeof(buf));
	p = buf;

	while (isspace(p[0])) {
		p++;
	}

	if (!p[0] || strncasecmp(p, "cd", 2) != 0 || (p[2] != '\0' && !isspace(p[2]))) {
		return -1;
	}

	if (read_fc511_block(in, data->NAD, data->NMD, data->ad, data->md)) {
		return -1;
	}

   	in.getline(buf, sizeof(buf));
	p = buf;

	while (isspace(p[0])) {
		p++;
	}

	if (!p[0] || strncasecmp(p, "cm", 2) != 0 || (p[2] != '\0' && !isspace(p[2]))) {
		return -1;
	}

	if (read_fc511_block(in, data->NAM, data->NMM, data->am, data->mm)) {
		return -1;
	}

	/* FIXME: maybe this is not the best place */
	c81_data_do_stall(data, dcltol);

   	return 0;
}

extern "C" int
c81_data_nrel_read(std::istream& in, c81_data* data, const doublereal dcltol)
{
   	char buf[128];	// 81 should suffice; let's make it 128

   	// header
   	in.getline(buf, sizeof(buf));

	// keep the first one
	memcpy(data->header, buf, sizeof(data->header));
	data->header[STRLENOF(data->header)] = '\0';

	// discard the second one
	in.getline(buf, sizeof(buf));

	int n;
	in >> n;
	if (n != 1) {
		silent_cerr("NREL format: can only handle one airfoil per file (got " << n << ")" << std::endl);
		return -1;
	}

	// discard remaining lines
	for (int i = 0; i < 12; i++) {
		// discard
		in.getline(buf, sizeof(buf));
	}

	std::streampos pos = in.tellg();
	for (n = 0; in; n++) {
		in.getline(buf, sizeof(buf));
	}

	in.clear();
	in.seekg(pos);
	n--;

	if (n < 2) {
		silent_cerr("insufficient number of data points n=" << n << " in NREL data format (n >= 2 required)" << std::endl);
		return -1;
	}

	/* lift */
	data->NML = 1;
	data->ml = new doublereal[data->NML];
	data->ml[0] = 0.;
	data->NAL = n;
	data->al = new doublereal[2*data->NAL];

	/* drag */
	data->NMD = 1;
	data->md = new doublereal[data->NMD];
	data->md[0] = 0.;
	data->NAD = n;
	data->ad = new doublereal[2*data->NAD];

	/* moment */
	data->NMM = 1;
	data->mm = new doublereal[data->NMM];
	data->mm[0] = 0.;
	data->NAM = n;
	data->am = new doublereal[2*data->NAM];

	for (int i = 0; i < n; i++) {
		doublereal alpha, cl, cd, cm;
		in >> alpha >> cl >> cd >> cm;

		if (i == 0 && alpha != -180.) {
			silent_cerr("warning: NREL data format: first alpha=" << alpha << " != -180.)" << std::endl);
		} else if (i == n - 1 && alpha != 180.) {
			silent_cerr("warning: NREL data format: last alpha=" << alpha << " != 180.)" << std::endl);
		}

		data->al[i] = alpha;
		data->al[data->NAL + i] = cl;
		data->ad[i] = alpha;
		data->ad[data->NAD + i] = cd;
		data->am[i] = alpha;
		data->am[data->NAM + i] = cm;
	}

	/* FIXME: maybe this is not the best place */
	c81_data_do_stall(data, dcltol);

   	return 0;
}

extern "C" int
c81_data_read_free_format(std::istream& in, c81_data* data, const doublereal dcltol)
{
   	char buf[128];	// 81 should suffice; let's make it 128

   	/* header */
   	in.getline(buf, sizeof(buf));
	if (strncasecmp(buf, "# FREE FORMAT", STRLENOF("# FREE FORMAT")) == 0) {
   		in.getline(buf, sizeof(buf));
	}

	char	*p = std::strchr(buf, ';');
	if (p == NULL) {
		return -1;
	}

	size_t	len = STRLENOF(data->header);
	if (size_t(p - buf) < len) {
		len = size_t(p - buf);
	}

   	strncpy(data->header, buf, len);
   	data->header[len] = '\0';

	std::istringstream s(++p);

   	s >> data->NML;
   	s >> data->NAL;
   	s >> data->NMD;
   	s >> data->NAD;
   	s >> data->NMM;
   	s >> data->NAM;

	if (data->NAM <= 0 || data->NMM <= 0
			|| data->NAD <= 0 || data->NMD <= 0
			|| data->NAL <= 0 || data->NML <= 0) {
		return -1;
	}

   	/* lift */
   	data->ml = new doublereal[data->NML];
	for (int c = 0; c < data->NML; c++) {
   		in >> data->ml[c];
	}
	if (check_vec(data->ml, data->NML)) {
		return -1;
	}

   	data->al = new doublereal[(data->NML + 1)*data->NAL];
	for (int r = 0; r < data->NAL; r++) {
		/* NOTE: "<=" because the number of columns is data->NML + 1
		 * for the angle of attack */
		for (int c = 0; c <= data->NML; c++) {
			in >> data->al[r + data->NAL*c];
		}
	}
	if (check_vec(data->al, data->NAL)) {
		return -1;
	}

   	/* drag */
   	data->md = new doublereal[data->NMD];
	for (int c = 0; c < data->NMD; c++) {
   		in >> data->md[c];
	}
	if (check_vec(data->md, data->NMD)) {
		return -1;
	}

   	data->ad = new doublereal[(data->NMD + 1)*data->NAD];
	for (int r = 0; r < data->NAD; r++) {
		for (int c = 0; c <= data->NMD; c++) {
			in >> data->ad[r + data->NAD*c];
		}
	}
	if (check_vec(data->ad, data->NAD)) {
		return -1;
	}

   	/* moment */
   	data->mm = new doublereal[data->NMM];
	for (int c = 0; c < data->NMM; c++) {
   		in >> data->mm[c];
	}
	if (check_vec(data->mm, data->NMM)) {
		return -1;
	}

   	data->am = new doublereal[(data->NMM + 1)*data->NAM];
	for (int r = 0; r < data->NAM; r++) {
		for (int c = 0; c <= data->NMM; c++) {
			in >> data->am[r + data->NAM*c];
		}
	}
	if (check_vec(data->am, data->NAM)) {
		return -1;
	}

	/* FIXME: maybe this is not the best place */
	c81_data_do_stall(data, dcltol);

   	return 0;
}

extern "C" int
c81_data_write_free_format(std::ostream& out, c81_data* data)
{
	if (data == 0) {
		return -1;
	}

	out << "# FREE FORMAT" << std::endl;

   	out << data->header << ";"
     		<< " " << data->NML
     		<< " "<< data->NAL
     		<< " "<< data->NMD
     		<< " "<< data->NAD
     		<< " "<< data->NMM
     		<< " "<< data->NAM
		<< std::endl;

	/* lift */
	for (int c = 0; c < data->NML; c++) {
		out << data->ml[c] << " ";
	}
	out << std::endl;
	for (int r = 0; r < data->NAL; r++) {
		/* NOTE: "<=" because the number of columns is data->NML + 1
		 * for the angle of attack */
		for (int c = 0; c <= data->NML; c++) {
			out << data->al[r + data->NAL*c] << " ";
		}
		out << std::endl;
	}

	/* drag */
	for (int c = 0; c < data->NMD; c++) {
		out << data->md[c] << " ";
	}
	out << std::endl;
	for (int r = 0; r < data->NAD; r++) {
		for (int c = 0; c <= data->NMD; c++) {
			out << data->ad[r + data->NAD*c] << " ";
		}
		out << std::endl;
	}

	/* moment */
	for (int c = 0; c < data->NMM; c++) {
		out << data->mm[c] << " ";
	}
	out << std::endl;
	for (int r = 0; r < data->NAM; r++) {
		for (int c = 0; c <= data->NMM; c++) {
			out << data->am[r + data->NAM*c] << " ";
		}
		out << std::endl;
	}

   	return 0;
}

/* data array */
static int c81_ndata = 0;
static c81_data** __c81_pdata = NULL;

extern "C" c81_data*
get_c81_data(long int jpro)
{
   	if (__c81_pdata == NULL) {
      		return NULL;
   	}
   	return __c81_pdata[jpro];
}

extern "C" int
c81_data_set(long int jpro, c81_data* data)
{
   	if (__c81_pdata == NULL || jpro >= c81_ndata) {
      		c81_data** pp = NULL;
      		int ndata = c81_ndata;

      		if (jpro >= c81_ndata) {
	 		ndata = jpro+1;
      		}

      		pp = new c81_data*[ndata];
      		if (__c81_pdata != NULL) {
	 		for (int i = 0; i < ndata; i++) {
	    			pp[i] = __c81_pdata[i];
	 		}
	 		delete[] __c81_pdata;
      		}
      		__c81_pdata = pp;
      		c81_ndata = ndata;
   	} else if (__c81_pdata[jpro] != NULL) {
      		silent_cerr("profile " << jpro << "already defined" << std::endl);
      		return -1;
   	}

   	__c81_pdata[jpro] = data;

   	return 0;
}

/*
 * sistema i dati di stallo
 */
static int
do_stall(int NM, int NA, doublereal *a, doublereal *stall, const doublereal dcltol)
{
	for (int nm = 0; nm < NM; nm++) {
		int start = NA*(nm + 1);
		int na = NA/2;	// mid point

		// look for zero
		while (a[na] > 0.) {
			na--;
			if (na <= 0) {
				silent_cerr("do_stall: "
					"unable to find negative alpha range "
					"for Mach #" << nm + 1 << std::endl);
				return -1;
			}
		}

		while (a[na] < 0.) {
			na++;
			if (na >= NA) {
				silent_cerr("do_stall: "
					"unable to find positive alpha range "
					"for Mach #" << nm + 1 << std::endl);
				return -1;
			}
		}

		doublereal a0 = a[na];
		doublereal dcl0 = a[start + na];

		doublereal dcla0 = (a[start + na + 1] - dcl0)/(a[na + 1] - a0);

		/* cerca il punto superiore in cui cessa la linearita' */
		stall[nm] = 1.;
		stall[NM + nm] = 1.;
		stall[2*NM + nm] = 0.;

		for (int i = na + 2; i < NA; i++) {
			doublereal dcla = (a[start + i] - dcl0)/(a[i] - a0);
			if (dcla - dcla0 < -dcla0*dcltol) {

				/* alpha di stallo superiore */
				stall[nm] = a[i - 1];

				/* mette temporaneamente il Cp di stallo */
				stall[2*NM + nm] = a[start + i - 1];
				break;
			}
		}

		dcla0 = (a[start + na - 1] - dcl0)/(a[na - 1] - a0);

		/* cerca il punto inferiore in cui cessa la linearita' */
		for (int i = na - 2; i >= 0; i--) {
			doublereal dcla = (a[start + i] - dcl0)/(a[i] - a0);
			if (dcla - dcla0 < -dcla0*dcltol) {

				/* stallo inferiore */
				stall[NM + nm] = a[i + 1];

				/* sottrae il cl allo stallo inferiore */
				stall[2*NM + nm] -= a[start + i + 1];

				/* divide per il delta alpha */
				stall[2*NM + nm] /= (stall[nm] - stall[NM + nm]);
				break;
			}
		}
	}

	return 0;
}

int
c81_data_do_stall(c81_data *data, const doublereal dcltol)
{
	if (data == NULL || data->NML <= 0 || data->NMM <= 0) {
		return -1;
	}

	data->stall = new doublereal[3*data->NML];
	do_stall(data->NML, data->NAL, data->al, data->stall, dcltol);

	data->mstall = new doublereal[3*data->NMM];
	do_stall(data->NMM, data->NAM, data->am, data->mstall, dcltol);

	return 0;
}

static int
flip_one(int NM, int NA, doublereal *v, int ss)
{
	int s = -1;
	for (int m = 0; m <= NM; m++) {
		if (m > 0) {
			s = ss;
		}

		for (int a = 0; a < NA/2; a++) {
			doublereal d = v[a];
			v[a] = s*v[NA - a - 1];
			v[NA - a - 1] = s*d;
		}

		// NA odd: change sign
		if (NA%2) {
			v[NA/2] *= s;
		}

		v += NA;
	}

	return 0;
}

extern "C" int 
c81_data_flip(c81_data *data)
{
	int rc;

	rc = flip_one(data->NML, data->NAL, data->al, -1);
	if (rc) {
		return rc;
	}

	rc = flip_one(data->NMD, data->NAD, data->ad, 1);
	if (rc) {
		return rc;
	}

	rc = flip_one(data->NMM, data->NAM, data->am, -1);
	return rc;
}

