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

int
do_grid(struct n2m_buffer *b, FILE **f)
{
	int ID, CP = 0;
	double X1, X2, X3;
	char form[1024];
	
	ID = get_int(b, NASTRAN_SECOND, NULL);
	CP = get_int(b, NASTRAN_THIRD, &CP);

	X1 = get_double(b, NASTRAN_FOURTH, NULL);
        X2 = get_double(b, NASTRAN_FIFTH, NULL);
        X3 = get_double(b, NASTRAN_SIXTH, NULL);
	
	fprintf(f[NASTRAN_FILE_OUT_REF],
		"reference: %8d, # GRID %d\n", reference_offset+ID, ID);
	if (CP == 0) {
		snprintf(form, sizeof(form), "\treference,   global");
	} else {
		snprintf(form, sizeof(form), "\treference, %8d", CP);
	}
	fprintf(f[NASTRAN_FILE_OUT_REF],
		"%s, %14.7e, %14.7e, %14.7e,\n", form, X1, X2, X3);
	fprintf(f[NASTRAN_FILE_OUT_REF], "%s, eye,\n", form);
	fprintf(f[NASTRAN_FILE_OUT_REF], "%s, null,\n", form);
	fprintf(f[NASTRAN_FILE_OUT_REF], "%s, null;\n", form);

	snprintf(form, sizeof(form), "\treference, %8d", ID);
	fprintf(f[NASTRAN_FILE_OUT_NODE],
		"structural: %8d, dynamic, # GRID %d\n",
		reference_offset+ID, ID);
	fprintf(f[NASTRAN_FILE_OUT_NODE], "%s, null,\n", form);
	fprintf(f[NASTRAN_FILE_OUT_NODE], "%s, eye,\n", form);
	fprintf(f[NASTRAN_FILE_OUT_NODE], "%s, null,\n", form);
	fprintf(f[NASTRAN_FILE_OUT_NODE], "%s, null;\n", form);

	return 0;
}

