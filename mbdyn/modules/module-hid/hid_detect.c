/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#include "stdlib.h"
#include "unistd.h"
#include "stdio.h"
#include "errno.h"
#include "stdint.h"
#include "string.h"

#define LO(x) ((x) & 0x0F)
#define HI(x) (((x) & 0xF0) >> 8)
#define HEX(x) ("0123456789abcdef"[(x)])

static const char *
type2str(uint8_t type)
{
	switch (type) {
	case 1:
		return "button";
	case 2:
		return "linctl";
	}
	return "unknwn";
}

int
main(int argc, char *argv[])
{
	FILE *fd;

	const char *device = "/dev/input/js0";
	
	if (argc > 1) {
		device = argv[1];
	}

	fd = fopen(device, "r");
	if (fd == NULL) {
		int save_errno = errno;
		fprintf(stderr, "fopen(\"%s\")=%d %s\n", device, save_errno, strerror(save_errno));
		exit(EXIT_FAILURE);
	}

	int nbt = 0, nlc = 0;
	uint8_t maxbt = 0, maxlc = 0;
	int preamble = 1;

	for (;;) {
		char buf[8];
		uint32_t cnt;
		int16_t value;
		uint8_t type, idx;

		size_t rc;

		rc = fread((void *)&buf[0], 1, sizeof(buf), fd);
		if (rc < 1) {
			fprintf(stderr, "fread returned no data\n");
			break;
		}

		cnt = *((uint32_t *)&buf[0]);
		value = *((int16_t *)&buf[4]);
		type = *((uint8_t *)&buf[6]);
		idx = *((uint8_t *)&buf[7]);

		if (preamble) {
			if (!(type & 0x80U)) {
				printf("buttons=%d (max=%u) linear controls=%d (max=%u)\n\n", nbt, maxbt, nlc, maxlc);
				preamble = 0;

			} else {
				switch (type & 0x7FU) {
				case 1:
					nbt++;
					if (idx > maxbt) {
						maxbt = idx;
					}
					break;

				case 2:
					nlc++;
					if (idx > maxlc) {
						maxlc = idx;
					}
					break;
				}

				printf("0x%8x %s[%2u]=%6d\n", cnt, type2str(type & 0x7FU), idx, value);

				continue;
			}
		}

		printf("0x%8x %s[%2u]=%6d\n", cnt, type2str(type), idx, value);
	}

	return 0;
}
