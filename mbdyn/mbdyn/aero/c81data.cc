/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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

#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <ac/f2c.h>
#include <cmath>

extern "C" {
#include <aerodc81.h>

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
do_c81_data_stall(c81_data *data, const doublereal dcltol);
static int
do_c81_stall(int NM, int NA, doublereal *a, doublereal *stall, const doublereal dcltol);

int
read_c81_data_free_format(std::istream& in, c81_data* data, const doublereal dcltol);

static int
get_int(const char *const from, int &i)
{
#ifdef HAVE_STRTOL
	char *endptr = NULL;
	i = strtol(from, &endptr, 10);
	if (endptr != NULL && endptr[0] != '\0' && !isspace(endptr[0])) {
		return -1;
	}
#else /* !HAVE_STRTOL */
   	i = atoi(buf);
#endif /* !HAVE_STRTOL */
	return 0;
}

static int
get_long(const char *const from, long &l)
{
#ifdef HAVE_STRTOL
	char *endptr = NULL;
	l = strtol(from, &endptr, 10);
	if (endptr != NULL && endptr[0] != '\0' && !isspace(endptr[0])) {
		return -1;
	}
#else /* !HAVE_STRTOL */
   	l = atol(buf);
#endif /* !HAVE_STRTOL */
	return 0;
}

static int
get_double(const char *const from, doublereal &d)
{
#ifdef HAVE_STRTOD
	char *endptr = NULL;
	d = strtod(from, &endptr);
	if (endptr != NULL && endptr[0] != '\0' && !isspace(endptr[0])) {
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

void
destroy_c81_data(c81_data* data)
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

int
read_c81_data(std::istream& in, c81_data* data, const doublereal dcltol)
{
   	char buf[BUFSIZ];	// 81 should suffice
   
   	/* header */
   	in.getline(buf, sizeof(buf));
	if (strncasecmp(buf, "# FREE FORMAT", STRLENOF("# FREE FORMAT")) == 0) {
		return read_c81_data_free_format(in, data, dcltol);
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
   
   	data->al = new doublereal[(data->NML+1)*data->NAL];
   	get_c81_mat(in, data->al, data->NAL, data->NML+1);
   
   	/* drag */
   	data->md = new doublereal[data->NMD];
   	get_c81_vec(in, data->md, data->NMD);
   
   	data->ad = new doublereal[(data->NMD+1)*data->NAD];      
   	get_c81_mat(in, data->ad, data->NAD, data->NMD+1);
   
   	/* moment */
   	data->mm = new doublereal[data->NMM];
   	get_c81_vec(in, data->mm, data->NMM);
   
   	data->am = new doublereal[(data->NMM+1)*data->NAM];
   	get_c81_mat(in, data->am, data->NAM, data->NMM+1);

	/* FIXME: maybe this is not the best place */
	do_c81_data_stall(data, dcltol);
   
   	return 0;
}

int
write_c81_data(std::ostream& out, c81_data* data)
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

int
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

int
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

int
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

int
read_fc511_data(std::istream& in, c81_data* data, const doublereal dcltol)
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
	do_c81_data_stall(data, dcltol);
   
   	return 0;
}

int
read_c81_data_free_format(std::istream& in, c81_data* data, const doublereal dcltol)
{
   	char buf[128];	// 81 should suffice; let's make it 128
   
   	/* header */
   	in.getline(buf, sizeof(buf));
	if (strncasecmp(buf, "# FREE FORMAT", STRLENOF("# FREE FORMAT")) == 0) {
   		in.getline(buf, sizeof(buf));
	}

	char	*p = strchr(buf, ';');
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
   
   	data->al = new doublereal[(data->NML + 1)*data->NAL];
	for (int r = 0; r < data->NAL; r++) {
		/* NOTE: "<=" because the number of columns is data->NML + 1
		 * for the angle of attack */
		for (int c = 0; c <= data->NML; c++) {
			in >> data->al[r + data->NAL*c];
		}
	}
   
   	/* drag */
   	data->md = new doublereal[data->NMD];
	for (int c = 0; c < data->NMD; c++) {
   		in >> data->md[c];
	}
   
   	data->ad = new doublereal[(data->NMD + 1)*data->NAD];      
	for (int r = 0; r < data->NAD; r++) {
		for (int c = 0; c <= data->NMD; c++) {
			in >> data->ad[r + data->NAD*c];
		}
	}
   
   	/* moment */
   	data->mm = new doublereal[data->NMM];
	for (int c = 0; c < data->NMM; c++) {
   		in >> data->mm[c];
	}
   
   	data->am = new doublereal[(data->NMM + 1)*data->NAM];
	for (int r = 0; r < data->NAM; r++) {
		for (int c = 0; c <= data->NMM; c++) {
			in >> data->am[r + data->NAM*c];
		}
	}

	/* FIXME: maybe this is not the best place */
	do_c81_data_stall(data, dcltol);
   
   	return 0;
}

int
write_c81_data_free_format(std::ostream& out, c81_data* data)
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

c81_data*
get_c81_data(long int jpro)
{
   	if (__c81_pdata == NULL) {
      		return NULL;
   	}
   	return __c81_pdata[jpro];
}

int 
set_c81_data(long int jpro, c81_data* data)
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
do_c81_stall(int NM, int NA, doublereal *a, doublereal *stall, const doublereal dcltol)
{
	for (int nm = 0; nm < NM; nm++) {
		int start = NA*(nm + 1);
		int na = NA/2;	// mid point

		// look for zero
		while (a[na] > 0.) {
			na--;
			if (na <= 0) {
				silent_cerr("do_c81_stall: "
					"unable to find negative alpha range "
					"for Mach #" << nm + 1 << std::endl);
				return -1;
			}
		}

		while (a[na] < 0.) {
			na++;
			if (na >= NA) {
				silent_cerr("do_c81_stall: "
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

static int
do_c81_data_stall(c81_data *data, const doublereal dcltol)
{
	if (data == NULL || data->NML <= 0 || data->NMM <= 0) {
		return -1;
	}

	data->stall = new doublereal[3*data->NML];
	do_c81_stall(data->NML, data->NAL, data->al, data->stall, dcltol);

	data->mstall = new doublereal[3*data->NMM];
	do_c81_stall(data->NMM, data->NAM, data->am, data->mstall, dcltol);

	return 0;
}

} // extern "C"
