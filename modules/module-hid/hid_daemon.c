/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
#include "time.h"
#include "pthread.h"

#define LO(x) ((x) & 0x0F)
#define HI(x) (((x) & 0xF0) >> 8)
#define HEX(x) ("0123456789abcdef"[(x)])
#define IS_BUTTON(type) (((type) & 0x7FU) == 1)
#define IS_LC(type) (((type) & 0x7FU) == 2)

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

typedef struct data_t {
	uint32_t cnt;
	int16_t value;
#define BUTTON_CLEARED 0x0U
#define BUTTON_SET 0x1U
#define BUTTON_CLEAR 0x2U
	uint8_t type;
	uint8_t idx;
} data_t;

typedef struct ptarg_t {
	pthread_t thread;
	pthread_mutex_t mutex;
	FILE *fd;

	int nbt;
	int nlc;
	int ndata;
	data_t *data;
} ptarg_t;

static int
read_one(FILE *fd, uint32_t *cnt, int16_t *value, uint8_t *type, uint8_t *idx)
{
	char buf[8];
	size_t rc;

	rc = fread((void *)&buf[0], 1, sizeof(buf), fd);
	if (rc < sizeof(buf)) {
		fprintf(stderr, "fread returned no data\n");
		return 1;
	}

	*cnt = *((uint32_t *)&buf[0]);
	*value = *((int16_t *)&buf[4]);
	*type = *((uint8_t *)&buf[6]);
	*idx = *((uint8_t *)&buf[7]);

	return 0;
}

static void *
purge_hid(void *arg)
{
	ptarg_t *ptarg = (ptarg_t *)arg;

	for (;;) {
		uint32_t cnt;
		int16_t value;
		uint8_t type, idx;
		int rc;

		rc = read_one(ptarg->fd, &cnt, &value, &type, &idx);
		if (rc) {
			/* error */
		}

		rc = idx;
		pthread_mutex_lock(&ptarg->mutex);
		if (IS_BUTTON(type)) {
			rc += ptarg->nlc;
			if (value) {
				ptarg->data[rc].value = BUTTON_SET;
			} else {
				ptarg->data[rc].value |= BUTTON_CLEAR;
			}

		} else {
			ptarg->data[rc].value = value;
		}
		ptarg->data[rc].cnt = cnt;
		ptarg->data[rc].type = type;
		ptarg->data[rc].idx = idx;
		pthread_mutex_unlock(&ptarg->mutex);
	}

	return NULL;
}

int
main(int argc, char *argv[])
{

	const char *device = "/dev/input/js0";
	long int dt = -1;
	int ndata = -1;

	struct timespec t0, t;
	ptarg_t ptarg = { 0 };
	uint8_t maxbt = 0, maxlc = 0;
	int preamble = 1;

	while (1) {
		int opt = getopt(argc, argv, "N:T:");
		if (opt == -1) {
			break;
		}

		switch (opt) {
		case 'N': {
			char *next = 0;
			long l;
			l = strtol(optarg, &next, 10);
			if (next == optarg || l <= 0) {
				/* error */
				exit(EXIT_FAILURE);
			}
			ndata = l;
			} break;

		case 'T': {
			char *next = 0;
			dt = strtol(optarg, &next, 10);
			if (next == optarg || dt <= 0) {
				/* error */
				exit(EXIT_FAILURE);
			}
			} break;

		default:
			exit(EXIT_FAILURE);
		}
	}

	if (dt <= 0) {
		fprintf(stderr, "dt must be defined and positive\n");
		exit(EXIT_FAILURE);
	}

	if (optind == argc) {
		/* error */

	} else if (optind == argc - 1) {
		device = argv[optind];

	} else {
		/* error */
	}
	
	ptarg.fd = fopen(device, "r");
	if (ptarg.fd == NULL) {
		int save_errno = errno;
		fprintf(stderr, "fopen(\"%s\")=%d %s\n", device, save_errno, strerror(save_errno));
		exit(EXIT_FAILURE);
	}

	for (;;) {
		int c;

		if (preamble) {
			uint32_t cnt;
			int16_t value;
			uint8_t type, idx;

			if (ndata != 0) {
				c = read_one(ptarg.fd, &cnt, &value, &type, &idx);
			}

			if (ndata == 0 || !(type & 0x80U)) {
				printf("buttons=%d (max=%u) linear controls=%d (max=%u)\n\n", ptarg.nbt, maxbt, ptarg.nlc, maxlc);
				preamble = 0;

				clock_gettime(CLOCK_MONOTONIC, &t0);
				t = t0;

				pthread_mutex_init(&ptarg.mutex, NULL);
				ptarg.ndata = ptarg.nbt + ptarg.nlc;
				ptarg.data = (data_t *)malloc(ptarg.ndata*sizeof(data_t));
				for (c = 0; c < ptarg.nlc; c++) {
					ptarg.data[c].idx = c;
					ptarg.data[c].value = 0;
					ptarg.data[c].type = 2;
				}

				for (c = 0; c < ptarg.nbt; c++) {
					ptarg.data[ptarg.nlc + c].idx = c;
					ptarg.data[ptarg.nlc + c].value = 0;
					ptarg.data[ptarg.nlc + c].type = 1;
				}

				c = pthread_create(&ptarg.thread, NULL, purge_hid, &ptarg);

			} else {
				switch (type & 0x7FU) {
				case 1:
					ptarg.nbt++;
					if (idx > maxbt) {
						maxbt = idx;
					}
					break;

				case 2:
					ptarg.nlc++;
					if (idx > maxlc) {
						maxlc = idx;
					}
					break;
				}

				printf("0x%8x %s[%2u]=%6d\n", cnt, type2str(type & 0x7FU), idx, value);

				if (ndata > 0) {
					ndata--;
				}

				continue;
			}
		}

		t.tv_nsec += dt;
		if (t.tv_nsec >= 1000000000L) {
			t.tv_nsec -= 1000000000L;
			t.tv_sec++;
		}

		c = clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &t, NULL);

		pthread_mutex_lock(&ptarg.mutex);
		printf("time=%ld.%09ld\n", t.tv_sec, t.tv_nsec);
		for (c = 0; c < ptarg.ndata; c++) {
			if (IS_BUTTON(ptarg.data[c].type)) {
				int8_t set = ptarg.data[c].value & BUTTON_SET;
				int8_t clear = ptarg.data[c].value & BUTTON_CLEAR;
				printf("0x%08x %s[%2u]=%s%s\n",
					ptarg.data[c].cnt,
					type2str(ptarg.data[c].type),
					ptarg.data[c].idx,
					set ? "set" : "not set", clear ? ",to be cleared" : "");
				if (clear) {
					ptarg.data[c].value = BUTTON_CLEARED;
				}

			} else {
				printf("0x%08x %s[%2u]=%6d\n",
					ptarg.data[c].cnt,
					type2str(ptarg.data[c].type),
					ptarg.data[c].idx,
					ptarg.data[c].value);
			}
		}
		pthread_mutex_unlock(&ptarg.mutex);
	}

	return 0;
}
