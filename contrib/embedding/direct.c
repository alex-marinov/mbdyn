/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include "mbconfig.h"

#include "malloc.h"
#include "ac/f2c.h"
#include "invsolwrap.h"

int
main(void)
{
	const char *sIn = "ides_in";
	const char *sOut = "ides_out";

	unsigned uInFrcLabel = 97;
	unsigned uInKinLabel = 98;
	unsigned uOutLabel = 0;

	void *pS = mb_sol_create("direct.mbd", "output");

	if (mb_sol_prepare(pS) == 0) {
		doublereal *inKinbuf;
		doublereal *inFrcbuf;
		doublereal *outbuf;
		inKinbuf = malloc(sizeof(doublereal)*1);
		inFrcbuf = malloc(sizeof(doublereal)*1);
		outbuf = malloc(sizeof(doublereal)*2);
		mb_sol_setbufin(pS, uInKinLabel, 1, inKinbuf);
		mb_sol_setbufin(pS, uInFrcLabel, 1, inFrcbuf);
		mb_sol_setbufout(pS, uOutLabel, 2, outbuf);

		// write to I buffer
		inKinbuf[0] = 0.;
		inFrcbuf[0] = 0.;

		if (mb_sol_start(pS) == 0) {
			size_t t;
			// read from O buffer
			printf("%g %g\n", outbuf[0], outbuf[1]);

			for (t = 0; ; t++) {
				int b;

				// write to I buffer
				inKinbuf[0] = 0.;
				if (t == 100) {
					inFrcbuf[0] = 2000.;

				} else {
					inFrcbuf[0] = 0.;
				}

				b = mb_sol_advance(pS);

				// read from O buffer
				printf("%g %g\n", outbuf[0], outbuf[1]);

				if (b != 0) {
					break;
				}
			}
		}

		free(inKinbuf);
		free(inFrcbuf);
		free(outbuf);
	}

	mb_sol_destroy(pS);

	return 0;
}

