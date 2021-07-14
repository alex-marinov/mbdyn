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

#include <n2m.h>

enum {
	CORD2R_P1 = 0,
	CORD2R_P2,
	CORD2R_P3,
	CORD2R_Q1,
	CORD2R_Q2,
	CORD2R_Q3,
	CORD2R_R1,
	CORD2R_R2,
	CORD2R_R3
};

int
do_cord2r_cont(struct n2m_buffer *b, FILE **f, struct n2m_buffer *form, double *coords)
{
	double R1, R2, R3;
	R1 = get_double(b, NASTRAN_SECOND, NULL);
        R2 = get_double(b, NASTRAN_THIRD, NULL);
        R3 = get_double(b, NASTRAN_FOURTH, NULL);
	
	fprintf(f[NASTRAN_FILE_OUT_REF],
		"%s, 3, %14.7e, %14.7e, %14.7e,\n",
		form->buf, 
		coords[CORD2R_Q1]-coords[CORD2R_P1],
		coords[CORD2R_Q2]-coords[CORD2R_P2],
		coords[CORD2R_Q3]-coords[CORD2R_P3]);
	fprintf(f[NASTRAN_FILE_OUT_REF],
		"\t                   1, %14.7e, %14.7e, %14.7e,\n",
		R1-coords[CORD2R_P1],
		R2-coords[CORD2R_P2],
		R3-coords[CORD2R_P3]);
	fprintf(f[NASTRAN_FILE_OUT_REF], "%s, null,\n", form->buf);
	fprintf(f[NASTRAN_FILE_OUT_REF], "%s, null;\n", form->buf);
	
	return 0;
}

int
do_cord2r(struct n2m_buffer *b, FILE **f, struct n2m_buffer *form, double *coords)
{
	int ID, CID = 0;
	
	ID = get_int(b, NASTRAN_SECOND, NULL);
        CID = get_int(b, NASTRAN_THIRD, NULL);
		
	coords[CORD2R_P1] = get_double(b, NASTRAN_FOURTH, NULL);
	coords[CORD2R_P2] = get_double(b, NASTRAN_FIFTH, NULL);
	coords[CORD2R_P3] = get_double(b, NASTRAN_SIXTH, NULL);

	coords[CORD2R_Q1] = get_double(b, NASTRAN_SEVENTH, NULL);
	coords[CORD2R_Q2] = get_double(b, NASTRAN_EIGHTH, NULL);
	coords[CORD2R_Q3] = get_double(b, NASTRAN_NINTH, NULL);
	
	fprintf(f[NASTRAN_FILE_OUT_REF],
		"reference: %7d, # CORD2R %d\n", ID, ID);
	if (CID == 0) {
		snprintf(form->buf, sizeof(form->buf),
			"\treference, global");
	} else {
		snprintf(form->buf, sizeof(form->buf),
			"\treference, %7d", CID);
	}
	fprintf(f[NASTRAN_FILE_OUT_REF], "%s, %14.7e, %14.7e, %14.7e,\n",
		form->buf,
		coords[CORD2R_P1], coords[CORD2R_P2], coords[CORD2R_P3]);

	return 0;
}

