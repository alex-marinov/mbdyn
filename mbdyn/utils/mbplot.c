/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <termios.h>

/* haldler to piped gnuplot */
FILE *gnuplot = NULL;

/* Use this variable to remember original terminal attributes. */
int interactive = 0;
struct termios saved_attributes;

void
cleanup(void)
{
	if (gnuplot != NULL) {
		pclose(gnuplot);
	}

	if (interactive) {
		tcsetattr(STDIN_FILENO, TCSANOW, &saved_attributes);
	}
}

void
set_input_mode(void)
{
	struct termios tattr;

	/* Make sure stdin is a terminal. */
	if (!isatty(STDIN_FILENO)) {
		fprintf(stderr, "Not a terminal.\n");
		exit(EXIT_FAILURE);
	}

	/* Save the terminal attributes so we can restore them later. */
	tcgetattr(STDIN_FILENO, &saved_attributes);

	/* Set the funny terminal modes. */
	tcgetattr(STDIN_FILENO, &tattr);
	tattr.c_lflag &= ~(ICANON|ECHO); /* Clear ICANON and ECHO. */
	tattr.c_cc[VMIN] = 0;
	tattr.c_cc[VTIME] = 0;
	if (tcsetattr(STDIN_FILENO, TCSAFLUSH, &tattr) < 0) {
		perror("tcsetattr");
		exit(EXIT_FAILURE);
	}

	interactive = 1;
}


int
main(int argc, char *argv[])
{
	char *cmd = NULL;
	unsigned int sleeptime = 1;
	int skip = 0;

	setbuf(stdout, NULL);

	while (1) {
		int opt = getopt(argc, argv, "t:");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 't':
			sleeptime = atoi(optarg);
			break;
		}
	}

	if (optind == argc) {
		fprintf(stderr, "usage: mbplot [-t <sleeptime>] <command>\n");
		exit(EXIT_FAILURE);
	}

	cmd = argv[optind];

	fflush(NULL);
	gnuplot = popen("gnuplot", "w");
	if (gnuplot == NULL) {
		fprintf(stderr, "popen(\"gnuplot\") failed\n");
		exit(EXIT_FAILURE);
	}

	atexit(cleanup);
	set_input_mode();

	fprintf(gnuplot, "set data style l\n%s\n", cmd);
	fflush(gnuplot);
	for ( ; ; ) {
		sleep(sleeptime);

		if (interactive) {
			char c = '\0';

			switch (read(STDIN_FILENO, &c, 1)) {
			case 0:
				break;

			case -1:
				if (errno != EAGAIN) {
					fprintf(stderr, "input error\n");
				}
				break;

			case 1:
				fprintf(stdout, "\r                                                                                ");
				switch (c) {
				case 'p':
					fprintf(stdout, "\rpausing...");
					skip = 1;
					break;

				case 'c':
					fprintf(stdout, "\rcontinuing...");
					skip = 0;
					break;

				case 'k':
					fprintf(stdout, "\rbailing out...\n");
					exit(EXIT_SUCCESS);

				case 't': {
						  char buf[BUFSIZ] = {'\0'};
						  int cnt = 0;
						  int nc = 0;

						  while (1) {

							  nc = read(STDIN_FILENO, &buf[cnt], sizeof(buf)-cnt-1);
							  if (nc == -1 && errno != EAGAIN) {
								  break;

							  } else if (nc > 0) {
							  	cnt += nc;
							  	if (buf[cnt-1] == '\n') {
							  		buf[cnt-1] == '\0';
									break;
								}
							  }
						  }
  						  if (nc == -1 && errno != EAGAIN) {
							  break;
						  }

						  buf[sizeof(buf)-1] = '\0';
						  nc = atoi(buf);
						  if (nc > 0) {
							  sleeptime = nc;
  							  fprintf(stdout, "\rsetting sleeptime = %d", sleeptime);
						  }
						  break;
					  }
				}
				break;
			}
		}
		
		if (!skip) {
			fprintf(gnuplot, "replot\n");
			fflush(gnuplot);
		}
	}

	return EXIT_SUCCESS;
}

