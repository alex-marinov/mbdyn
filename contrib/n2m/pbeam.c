/******************************************************************************
 * 
 * NASTRAN 2 MBDyn bulk data conversion tool
 * 
 * Copyright (c) 2003
 * 
 * Author:         Pierangelo Masarati <masarati@aero.polimi.it>
 * Affiliation:    Dipartimento di Ingegneria Aerospaziale
 *                 Politecnico di Milano
 *                 via La Masa 34, 20156 Milano (Italy)
 *                 http://www.aero.polimi.it
 * 
 * This program is part of the MBDyn package
 * http://www.mbdyn.org
 * https://mbdyn.aero.polimi.it/~masarati/mbdyn-index.html
 * 
 *****************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <n2m.h>

static int
do_pbeam_section(int PID, int MID, struct n2m_buffer *b, FILE **f)
{
	double A, I1, I2, I12 = 0., J = 0., NSM = 0.;
	
	A = get_double(b, NASTRAN_FOURTH, NULL);
	
	I1 = get_double(b, NASTRAN_FIFTH, NULL);
	I2 = get_double(b, NASTRAN_SIXTH, NULL);
	I12 = get_double(b, NASTRAN_SEVENTH, &I12);
	J = get_double(b, NASTRAN_EIGHTH, &J);
	NSM = get_double(b, NASTRAN_NINTH, &NSM);
	
	fprintf(f[NASTRAN_FILE_OUT_ELEM],
		"\t#beam: /* PID = */ %8d, /* MID = */ %8d,\n",
		PID, MID);
	fprintf(f[NASTRAN_FILE_OUT_ELEM], "\tlinear elastic generic, diag,\n");
	fprintf(f[NASTRAN_FILE_OUT_ELEM],
		"\t%14.7e, %14.7e, %14.7e,\n", A, A, A);
	fprintf(f[NASTRAN_FILE_OUT_ELEM],
		"\t%14.7e, %14.7e, %14.7e,\n", J, I2, I1);
		
	return 0;
}

int
do_pbeam_cont(struct n2m_buffer *b, FILE **f, int *cont)
{
	(*cont)++;

	if (b->buf[0] != ' ') {
		*cont = 0;
		return 0;
	}

	switch (*cont) {
	case 2: {
		struct n2m_buffer SO;
		get_string(b, NASTRAN_SECOND, &SO, "");
		if (strncasecmp(SO.buf, "YES", 3) == 0
			|| strcasecmp(SO.buf, "NO") == 0) {
			do_pbeam_section(-1, -1, b, f);
		}
		break;
	}
	}
	
	return 0;
}

int
do_pbeam(struct n2m_buffer *b, FILE **f)
{
	int PID, MID;
	
	PID = get_int(b, NASTRAN_SECOND, NULL);
        MID = get_int(b, NASTRAN_THIRD, NULL);
		
	do_pbeam_section(PID, MID, b, f);
	
	return 0;
}

