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

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif /* HAVE_GETOPT_H */

#include <n2m.h>

int make_rigid_bodies = 0;
int reference_offset = 0;

int
main(int argc, char *const argv[])
{
	FILE *f[NASTRAN_FILE_LAST];
	struct n2m_buffer b;
	struct n2m_card_data cards[NASTRAN_CARD_LAST];
	int i;
	
	struct n2m_buffer form;		/* per continuaz. */
	double coords[6];		/* per cord2r */

	int type = -1;			/* tipo di card */
	int cont = 0;			/* n. riga continuazione */

	/* Initializes file handlers */
	f[NASTRAN_FILE_IN] = stdin;
	f[NASTRAN_FILE_OUT_CTL] = stdout;
	f[NASTRAN_FILE_OUT_REF] = stdout;
	f[NASTRAN_FILE_OUT_NODE] = stdout;
	f[NASTRAN_FILE_OUT_ELEM] = stdout;
	f[NASTRAN_FILE_OUT_ERR] = stderr;

	for (i = 0; i < NASTRAN_CARD_LAST; i++) {
		cards[i].num = 0;
	}

	/* Command-line parsing */
#ifdef HAVE_GETOPT
	while (1) {
		char optstring[] = "bc:e:i:n:o:r:R:";
		int opt;

		opt = getopt(argc, argv, optstring);

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'b':
			make_rigid_bodies++;
			break;
		case 'c':
			f[NASTRAN_FILE_OUT_CTL] = fopen(optarg, "w");
			if (f[NASTRAN_FILE_OUT_CTL] == NULL) {
				fprintf(f[NASTRAN_FILE_OUT_ERR],
					"### unable to open file '%s'\n",
					optarg);
				exit(EXIT_FAILURE);
			}
			break;
		case 'e':
			f[NASTRAN_FILE_OUT_ELEM] = fopen(optarg, "w");
			if (f[NASTRAN_FILE_OUT_ELEM] == NULL) {
				fprintf(f[NASTRAN_FILE_OUT_ERR],
					"### unable to open file '%s'\n",
					optarg);
				exit(EXIT_FAILURE);
			}
			break;
		case 'i':
			f[NASTRAN_FILE_IN] = fopen(optarg, "r");
			if (f[NASTRAN_FILE_IN] == NULL) {
				fprintf(f[NASTRAN_FILE_OUT_ERR],
					"### unable to open file '%s'\n",
					optarg);
				exit(EXIT_FAILURE);
			}
			break;
		case 'n':
			f[NASTRAN_FILE_OUT_NODE] = fopen(optarg, "w");
			if (f[NASTRAN_FILE_OUT_NODE] == NULL) {
				fprintf(f[NASTRAN_FILE_OUT_ERR],
					"### unable to open file '%s'\n",
					optarg);
				exit(EXIT_FAILURE);
			}
			break;
		case 'o': {
			char buf[1024];
			int l;

			l = strlen(optarg);
			strncpy(buf, optarg, sizeof(buf));

			strcpy(buf+l, ".ctl");
			f[NASTRAN_FILE_OUT_CTL] = fopen(buf, "w");
			if (f[NASTRAN_FILE_OUT_CTL] == NULL) {
				fprintf(f[NASTRAN_FILE_OUT_ERR],
					"### unable to open file '%s'\n", buf);
				exit(EXIT_FAILURE);
			}

			strcpy(buf+l, ".elm");
			f[NASTRAN_FILE_OUT_ELEM] = fopen(buf, "w");
			if (f[NASTRAN_FILE_OUT_ELEM] == NULL) {
				fprintf(f[NASTRAN_FILE_OUT_ERR],
					"### unable to open file '%s'\n", buf);
				exit(EXIT_FAILURE);
			}
			
			strcpy(buf+l, ".nod");
			f[NASTRAN_FILE_OUT_NODE] = fopen(buf, "w");
			if (f[NASTRAN_FILE_OUT_NODE] == NULL) {
				fprintf(f[NASTRAN_FILE_OUT_ERR],
					"### unable to open file '%s'\n", buf);
				exit(EXIT_FAILURE);
			}

			strcpy(buf+l, ".ref");
			f[NASTRAN_FILE_OUT_REF] = fopen(buf, "w");
			if (f[NASTRAN_FILE_OUT_REF] == NULL) {
				fprintf(f[NASTRAN_FILE_OUT_ERR],
					"### unable to open file '%s'\n", buf);
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'r':
			f[NASTRAN_FILE_OUT_REF] = fopen(optarg, "w");
			if (f[NASTRAN_FILE_OUT_REF] == NULL) {
				fprintf(f[NASTRAN_FILE_OUT_ERR],
					"### unable to open file '%s'\n",
					optarg);
				exit(EXIT_FAILURE);
			}
			break;
		case 'R':
			reference_offset = atoi(optarg);
			break;
		}
	}
	argc -= optind-1;
	argv += optind-1;
#endif /* HAVE_GETOPT */

	/* Extra command-line parsing */
	if (argc > 1) {
		if (f[NASTRAN_FILE_IN] != stdin) {
			fprintf(f[NASTRAN_FILE_OUT_ERR],
				"### ignoring arg `%s'\n", argv[1]);
		} else {
			f[NASTRAN_FILE_IN] = fopen(argv[1], "r");
			if (f[NASTRAN_FILE_IN] == NULL) {
				fprintf(f[NASTRAN_FILE_OUT_ERR],
					"### unable to open file '%s'\n",
					argv[1]);
				exit(EXIT_FAILURE);
			}
		}
	}

	/* Main cycle */
	while (get_buf(f[NASTRAN_FILE_IN], &b)) {
		if (cont > 0) {
			int def = 0;
			
			if (b.buf[0] != ' ') {
				def = 1;
			}
			// card continuation
			switch (type) {
			case NASTRAN_CARD_CONM2:
				do_conm2_cont(&b, f, &form, def);
				cont = 0;
				break;
			case NASTRAN_CARD_CORD2R:
				do_cord2r_cont(&b, f, &form, coords);
				cont = 0;
				break;
			case NASTRAN_CARD_PBEAM:
				do_pbeam_cont(&b, f, &cont);
				break;
			default:
				fprintf(f[NASTRAN_FILE_OUT_ERR],
					"### unable to handle continuation card\n%s\nbailing out ...\n", 
					b.buf);
				exit(EXIT_FAILURE);
			}

			continue;
		}

		if (strncasecmp(b.buf, "CONM2", 5) == 0) {
			type = NASTRAN_CARD_CONM2;
			do_conm2(&b, f, &form);
			cont = 1;
		} else if (strncasecmp(b.buf, "CORD2R", 6) == 0) {
			type = NASTRAN_CARD_CORD2R;
			do_cord2r(&b, f, &form, coords);
			cont = 1;
		} else if (strncasecmp(b.buf, "GRID", 4) == 0) {
			type = NASTRAN_CARD_GRID;
			do_grid(&b, f);
		} else if (strncasecmp(b.buf, "PBEAM", 5) == 0) {
			type = NASTRAN_CARD_PBEAM;
			do_pbeam(&b, f);
			cont = 1;
		}

		cards[type].num++;
	}

	fprintf(f[NASTRAN_FILE_OUT_CTL],
		"\t#references: \t\t%8d; # CORD2R\n",
		cards[NASTRAN_CARD_CORD2R].num);
	fprintf(f[NASTRAN_FILE_OUT_CTL],
		"\tstructural nodes: \t%8d; # GRID\n",
		cards[NASTRAN_CARD_GRID].num);
	fprintf(f[NASTRAN_FILE_OUT_CTL],
		"\t%srigid bodies: \t\t%8d; # CONM2\n",
		(make_rigid_bodies ? "" : "#"),
		cards[NASTRAN_CARD_CONM2].num);

	return 0;
}

