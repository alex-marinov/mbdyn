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
	char form[1024];

	int cont = 0;
	
	if (argc > 1) {
		fin = fopen(argv[1], "r");
		if (fin == NULL) {
			fprintf(stderr, "### unable to open file '%s'\n", argv[1]);
			exit(EXIT_FAILURE);
		}
	}

	while (fgets(buf, sizeof(buf), fin)) {
		if (cont == 1) {
			if (buf[0] == ' ') {
				// card continuation
				char token[16];
				double I11 = 0., I21 = 0., I22 = 0., I31 = 0., I32 = 0., I33 = 0.;
				
				memset(token+8, 0, 8);
				strncpy(token, buf+8, 8);
				trim(token);
				double_exp(token);
				if (token[0] != '\0') {
					I11 = atof(token);
				}
				
				memset(token+8, 0, 8);
				strncpy(token, buf+16, 8);
				trim(token);
				double_exp(token);
				if (token[0] != '\0') {
					I21 = -atof(token);
				}
				
				memset(token+8, 0, 8);
				strncpy(token, buf+24, 8);
				trim(token);
				double_exp(token);
				if (token[0] != '\0') {
					I22 = atof(token);
				}
				
				memset(token+8, 0, 8);
				strncpy(token, buf+32, 8);
				trim(token);
				double_exp(token);
				if (token[0] != '\0') {
					I31 = -atof(token);
				}
				
				memset(token+8, 0, 8);
				strncpy(token, buf+40, 8);
				trim(token);
				double_exp(token);
				if (token[0] != '\0') {
					I32 = -atof(token);
				}
				
				memset(token+8, 0, 8);
				strncpy(token, buf+48, 8);
				trim(token);
				double_exp(token);
				if (token[0] != '\0') {
					I33 = atof(token);
				}
				
				fprintf(fout, "%s, diag, %e, %e, %e,\n", form, I11, I21, I31);
				fprintf(fout, "\t\t\t\t\t\t%e, %e,\n", I22, I32);
				fprintf(fout, "\t\t\t\t\t\t\t\t%e;\n", I33);
			} else {
				// no continuation
				fprintf(fout, "%s, null;\n", form);
			}
			
			cont = 0;
		}
		
		if (strncasecmp(buf, "CONM2", 5) == 0) {
			char token[16];
			int EID, G, CID;
			double M, X1, X2, X3;
			
			memset(token+8, 0, 8);
			strncpy(token, buf+8, 8);
			trim(token);
			EID = atoi(token); /* non viene usato */

			memset(token+8, 0, 8);
			strncpy(token, buf+16, 8);
			trim(token);
			G = atoi(token);

			memset(token+8, 0, 8);
			strncpy(token, buf+24, 8);
			trim(token);
			CID = atoi(token);

			memset(token+8, 0, 8);
			strncpy(token, buf+32, 8);
			trim(token);
			double_exp(token);
			M = atof(token);

			memset(token+8, 0, 8);
			strncpy(token, buf+40, 8);
			trim(token);
			double_exp(token);
			X1 = atof(token);

			memset(token+8, 0, 8);
			strncpy(token, buf+48, 8);
			trim(token);
			double_exp(token);
			X2 = atof(token);

			memset(token+8, 0, 8);
			strncpy(token, buf+56, 8);
			trim(token);
			double_exp(token);
			X3 = atof(token);
						
			fprintf(fout, "# EID=%d, G=%d\n", EID, G);
			fprintf(fout, "\t%e,\n", M);
			if (CID == 0) {
				snprintf(form, sizeof(form), "\treference, global");
			} else {
				snprintf(form, sizeof(form), "\treference, %d", CID);
			}
			fprintf(fout, "%s, %e, %e, %e,\n", form, X1, X2, X3);

			cont = 1;
		}

		memset(buf, 0, strlen(buf));
	}

	return 0;
}

