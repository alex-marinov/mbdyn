#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

extern int trim(char *s);
extern int double_exp(char *s);

int
main(int argc, const char *const argv[])
{
	FILE *fin = stdin;
	FILE *fout = stdout;
	char buf[1024];
	
	if (argc > 1) {
		fin = fopen(argv[1], "r");
		if (fin == NULL) {
			fprintf(stderr, "### unable to open file '%s'\n", argv[1]);
			exit(EXIT_FAILURE);
		}
	}

	while (fgets(buf, sizeof(buf), fin)) {
		if (strncasecmp(buf, "GRID", 4) == 0) {
			char token[16];
			int ID, CP;
			double X1, X2, X3;
			char form[1024];
			
			memset(token+8, 0, 8);
			strncpy(token, buf+8, 8);
			trim(token);
			ID = atoi(token);

			memset(token+8, 0, 8);
			strncpy(token, buf+16, 8);
			trim(token);
			if (token[0] == '\0') {
				CP = 0;
			} else {
				CP = atoi(token);
			}

			memset(token+8, 0, 8);
			strncpy(token, buf+24, 8);
			trim(token);
			double_exp(token);
			X1 = atof(token);

			memset(token+8, 0, 8);
			strncpy(token, buf+32, 8);
			trim(token);
			double_exp(token);
			X2 = atof(token);

			memset(token+8, 0, 8);
			strncpy(token, buf+40, 8);
			trim(token);
			double_exp(token);
			X3 = atof(token);
						
			fprintf(fout, "reference: %d,\n", ID);
			if (CP == 0) {
				snprintf(form, sizeof(form), "\treference, global");
			} else {
				snprintf(form, sizeof(form), "\treference, %d", CP);
			}
			fprintf(fout, "%s, %e, %e, %e,\n", form, X1, X2, X3);
			fprintf(fout, "%s, eye,\n", form);
			fprintf(fout, "%s, null,\n", form);
			fprintf(fout, "%s, null;\n", form);
		}
	}

	return 0;
}

