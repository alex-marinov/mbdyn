/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

#include <ac/iostream>
#include <ac/iomanip>
#include <ac/math.h>

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
do_c81_data_stall(c81_data *data);
static int
do_c81_stall(int NM, int NA, double *a, double *stall);

static int 
get_vec(std::istream& in, double* v, int nrows)
{
   	if (!in || v == NULL || nrows < 1) {
      		return -1;
   	}
   
   	for (int i = 0; i < nrows; i++) {
      		in >> v[i];
   	}
   
   	return 0;
}

static int 
get_mat(std::istream& in, double* m, int nrows, int ncols)
{
   	if (!in || m == NULL || nrows < 1 || ncols < 1) {
      		return -1;
   	}
   
   	for (int i = 0; i < nrows; i++) {
      		for (int j = 0; j < ncols; j++) {
	 		in >> m[i+nrows*j];
      		}
   	}
   
   	return 0;
}

static int
put_row(std::ostream& out, double* v, int dim, int ncols, int first = 0)
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
put_vec(std::ostream& out, double* v, int nrows)
{
   	if (!out || v == NULL || nrows < 1) {
      		return -1;
   	}

   	put_row(out, v, 1, nrows);
   
   	return 0;
}

static int
put_mat(std::ostream& out, double* m, int nrows, int ncols)
{
   	if (!out || m == NULL || nrows < 1 || ncols < 1) {
      		return -1;
   	}
     
   	for (int i = 0; i < nrows; i++) {
      		put_row(out, m+i, nrows, ncols, 1);
   	}
	 
   	return 0;
}

int
read_c81_data(std::istream& in, c81_data* data)
{
   	char buf[1024];      
   
   	/* header */
   	in.getline(buf, sizeof(buf));

#ifdef HAVE_STRTOL
   	buf[42] = '\0';

	char *endptr = NULL;
   	data->NAM = strtol(buf + 40, &endptr, 10);
	if (endptr != NULL && endptr[0] != '\0') {
		return -1;
	}
   	buf[40] = '\0';

   	data->NMM = strtol(buf + 38, &endptr, 10);
	if (endptr != NULL && endptr[0] != '\0') {
		return -1;
	}
   	buf[38] = '\0';
   
   	data->NAD = strtol(buf + 36, &endptr, 10);
	if (endptr != NULL && endptr[0] != '\0') {
		return -1;
	}
   	buf[36] = '\0';
   	data->NMD = strtol(buf + 34, &endptr, 10);
	if (endptr != NULL && endptr[0] != '\0') {
		return -1;
	}
   	buf[34] = '\0';
   
   	data->NAL = strtol(buf + 32, &endptr, 10);
	if (endptr != NULL && endptr[0] != '\0') {
		return -1;
	}
   	buf[32] = '\0';
   	data->NML = strtol(buf + 30, &endptr, 10);
	if (endptr != NULL && endptr[0] != '\0') {
		return -1;
	}
   	buf[30] = '\0';
#else /* !HAVE_STRTOL */
   	data->NAM = atoi(buf+40);
   	buf[40] = '\0';
   	data->NMM = atoi(buf+38);
   	buf[38] = '\0';
   
   	data->NAD = atoi(buf+36);
   	buf[36] = '\0';
   	data->NMD = atoi(buf+34);
   	buf[34] = '\0';
   
   	data->NAL = atoi(buf+32);
   	buf[32] = '\0';
   	data->NML = atoi(buf+30);
   	buf[30] = '\0';
#endif /* !HAVE_STRTOL */

	if (data->NAM <= 0 || data->NMM <= 0
			|| data->NAD <= 0 || data->NMD <= 0
			|| data->NAL <= 0 || data->NML <= 0) {
		return -1;
	}
   
   	strncpy(data->header, buf, 30);
   	data->header[30] = '\0';
   
   	/* lift */
   	data->ml = new double[data->NML];
   	get_vec(in, data->ml, data->NML);
   
   	data->al = new double[(data->NML+1)*data->NAL];
   	get_mat(in, data->al, data->NAL, data->NML+1);
   
   	/* drag */
   	data->md = new double[data->NMD];
   	get_vec(in, data->md, data->NMD);
   
   	data->ad = new double[(data->NMD+1)*data->NAD];      
   	get_mat(in, data->ad, data->NAD, data->NMD+1);
   
   	/* moment */
   	data->mm = new double[data->NMM];
   	get_vec(in, data->mm, data->NMM);
   
   	data->am = new double[(data->NMM+1)*data->NAM];
   	get_mat(in, data->am, data->NAM, data->NMM+1);

	/* FIXME: maybe this is not the best place */
	do_c81_data_stall(data);
   
   	return 0;
}

int
write_c81_data(std::ostream& out, c81_data* data)
{
   	out << data->header
     		<< std::setw(2) << data->NML
     		<< std::setw(2) << data->NAL
     		<< std::setw(2) << data->NMD
     		<< std::setw(2) << data->NAD
     		<< std::setw(2) << data->NMM
     		<< std::setw(2) << data->NAM
		<< std::endl;
   
   	put_vec(out, data->ml, data->NML);
   	put_mat(out, data->al, data->NAL, data->NML+1);
   
   	put_vec(out, data->md, data->NMD);
   	put_mat(out, data->ad, data->NAD, data->NMD+1);
   
   	put_vec(out, data->mm, data->NMM);
   	put_mat(out, data->am, data->NAM, data->NMM+1);
   
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
      		std::cerr << "profile " << jpro << "already defined" << std::endl;
      		return -1;
   	}
   
   	__c81_pdata[jpro] = data;
   
   	return 0;
}

/*
 * sistema i dati di stallo
 */
static const double dcptol = 1.e-2;

static int
do_c81_stall(int NM, int NA, double *a, double *stall)
{
	for (int nm = 0; nm < NM; nm++) {
		int start = NA*(nm+1);
		int na = NA/2;	/* punto di mezzo */
		double a0 = a[na];
		double dcp0 = a[start+na];
		double dcpa0 = (a[start+na+1]-dcp0)/(a[na+1]-a0);

		/* cerca il punto superiore in cui cessa la linearita' */
		stall[nm] = 1.;
		stall[NM+nm] = 1.;
		stall[2*NM+nm] = 0.;

		for (int i = na+2; i < NA; i++) {
			double dcpa = (a[start+i]-dcp0)/(a[i]-a0);
			if (fabs(dcpa-dcpa0) > dcptol) {

				/* alpha di stallo superiore */
				stall[nm] = a[i-1];

				/* mette temporaneamente il Cp di stallo */
				stall[2*NM+nm] = 
					a[start+i-1];
				break;
			}
		}

		dcpa0 = (a[start+na-1]-dcp0)/(a[na-1]-a0);
		
		/* cerca il punto inferiore in cui cessa la linearita' */
		for (int i = na-2; i >= 0; i--) {
			double dcpa = (a[start+i]-dcp0)/(a[i]-a0);
			if (fabs(dcpa-dcpa0) > dcptol) {
			
				/* stallo inferiore */
				stall[NM+nm] = a[i+1];

				/* sottrae il cp allo stallo inferiore */
				stall[2*NM+nm] -= 
					a[start+i+1];
					
				/* divide per il delta alpha */
				stall[2*NM+nm] /= (stall[nm] - stall[NM+nm]);
				break;
			}
		}
	}

	return 0;
}

static int
do_c81_data_stall(c81_data *data)
{
	if (data == NULL || data->NML <= 0 || data->NMM <= 0) {
		return -1;
	}

	data->stall = new double[3*data->NML];
	do_c81_stall(data->NML, data->NAL, data->al, data->stall);

	data->mstall = new double[3*data->NMM];
	do_c81_stall(data->NMM, data->NAM, data->am, data->mstall);

	return 0;
}

} /* extern "C" */
