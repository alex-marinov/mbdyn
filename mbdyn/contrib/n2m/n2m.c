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

int do_conm2(struct n2m_buffer *b, FILE *fout, struct n2m_buffer *form);
int do_conm2_cont(struct n2m_buffer *b, FILE *fout, struct n2m_buffer *form, int def);

int
main(int argc, const char *const argv[])
{
	FILE *fin = stdin;
	FILE *fout_node = stdout;
	FILE *fout_elem = stdout;
	struct n2m_buffer b;
	
	struct n2m_buffer form;		/* per continuaz. */
	double coords[6];		/* per cord2r */

	int type = -1;
	int cont = 0;
	
	if (argc > 1) {
		fin = fopen(argv[1], "r");
		if (fin == NULL) {
			fprintf(stderr,
				"### unable to open file '%s'\n", argv[1]);
			exit(EXIT_FAILURE);
		}
	}

	while (get_buf(fin, &b)) {
		if (cont > 0) {
			int def = 0;
			
			if (b.buf[0] != ' ') {
				def = 1;
			}
			// card continuation
			switch (type) {
			case NASTRAN_CARD_CONM2:
				do_conm2_cont(&b, fout_elem, &form, def);
				cont = 0;
				break;
			case NASTRAN_CARD_CORD2R:
				do_cord2r_cont(&b, fout_elem, &form, coords);
				cont = 0;
				break;
			default:
				fprintf(stderr,
					"### unable to handle continuation card\n%s\nbailing out ...\n", 
					b.buf);
				exit(EXIT_FAILURE);
			}

			continue;
		}

		if (strncasecmp(b.buf, "CONM2", 5) == 0) {
			type = NASTRAN_CARD_CONM2;
			do_conm2(&b, fout_elem, &form);
			cont = 1;
		} else if (strncasecmp(b.buf, "CORD2R", 6) == 0) {
			type = NASTRAN_CARD_CORD2R;
			do_cord2r(&b, fout_elem, &form, coords);
			cont = 1;
		} else if (strncasecmp(b.buf, "GRID", 4) == 0) {
			type = NASTRAN_CARD_GRID;
			do_grid(&b, fout_node);
		}
	}

	return 0;
}

