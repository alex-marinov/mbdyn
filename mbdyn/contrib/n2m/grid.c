/******************************************************************************
 * 
 * NASTRAN 2 MBDyn bulk data conversion tool
 * 
 * Copyright (c) 2000
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

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <n2m.h>

int
do_grid(struct n2m_buffer *b, FILE *fout)
{
	int ID, CP = 0;
	double X1, X2, X3;
	char form[1024];
	
	ID = get_int(b, NASTRAN_SECOND, NULL);
	CP = get_int(b, NASTRAN_THIRD, &CP);

	X1 = get_double(b, NASTRAN_FOURTH, NULL);
        X2 = get_double(b, NASTRAN_FIFTH, NULL);
        X3 = get_double(b, NASTRAN_SIXTH, NULL);
	
	fprintf(fout, "reference: %8d, # GRID %d\n", ID, ID);
	if (CP == 0) {
		snprintf(form, sizeof(form), "\treference,   global");
	} else {
		snprintf(form, sizeof(form), "\treference, %8d", CP);
	}
	fprintf(fout, "%s, %14.7e, %14.7e, %14.7e,\n", form, X1, X2, X3);
	fprintf(fout, "%s, eye,\n", form);
	fprintf(fout, "%s, null,\n", form);
	fprintf(fout, "%s, null;\n", form);

	return 0;
}

