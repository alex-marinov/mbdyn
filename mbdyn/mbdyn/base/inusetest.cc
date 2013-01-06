/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <stdlib.h>
#include <iostream>
#include <cstring>
#include "ac/getopt.h"

#ifdef USE_MULTITHREAD

#include <unistd.h>
#include "ac/pthread.h"		/* includes POSIX semaphores */

#include "mbconfig.h"
#include "veciter.h"
#include "filename.h"

unsigned dst = 100;
unsigned rst = 10;

class A : public InUse {
private:
	unsigned label;
	unsigned who;
public:
	A(unsigned l) : label(l) {};
	virtual ~A(void) {};

	unsigned Label(void) const { return label; };
	void Set(unsigned w) {
		who = w;

		unsigned s = dst;
		if (rst) {
			s += rand()%rst;
		}
		if (s) {
			usleep(s);
		}
	};
	unsigned Get(void) const { return who; };
};

struct Arg {
	unsigned n;
	unsigned cnt;
	bool stop;
	pthread_t t;
	sem_t s;
	MT_VecIter<A *> i;
	unsigned *c;
};

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;

static void
f2(VecIter<A *> &i, unsigned n, unsigned &cnt)
{
	A* pA = NULL;

	cnt = 0;

	if (i.bGetFirst(pA)) {
		do {
			pA->Set(n);
			cnt++;
		} while (i.bGetNext(pA));
	}
}

void *
f(void *p)
{
	Arg *arg = (Arg *)p;

	while (true) {
		sem_wait(&arg->s);

		if (arg->stop) {
			break;
		}

		f2(arg->i, arg->n, arg->cnt);

		pthread_mutex_lock(&mutex);
		--*(arg->c);
		silent_cout("f: count " << *arg->c + 1 << " => " << *arg->c << std::endl);
		if (*arg->c == 0) {
			silent_cout("f: finished" << std::endl);
			pthread_cond_signal(&cond);
		}
		pthread_mutex_unlock(&mutex);
	}

	return NULL;
}

int
main(int argc, char* argv[])
{
	unsigned size = 1000;
	unsigned nt = 1;
	unsigned loops = 1;

	if (argc == 1) {
usage:;
		char *s = std::strrchr(argv[0], DIR_SEP);

		if (s) {
			s++;
		} else {
			s = argv[0];
		}

		std::cout << "usage: " << s << " [lnsSt]" << std::endl
			<< "\t-l <loops>" << std::endl
			<< "\t-n <size>" << std::endl
			<< "\t-s <sleeptime>" << std::endl
			<< "\t-S <random sleeptime>" << std::endl
			<< "\t-t <threads number>" << std::endl;
		exit(EXIT_SUCCESS);
	}

	while (true) {
		char	*next;
		int	opt = getopt(argc, argv, "l:n:s:S:t:");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'l':
			loops = strtoul(optarg, &next, 10);
			break;

		case 'n':
			size = strtoul(optarg, &next, 10);
			break;

		case 's':
			dst = strtoul(optarg, &next, 10);
			break;

		case 'S':
			rst = strtoul(optarg, &next, 10);
			break;

		case 't':
			nt = strtoul(optarg, &next, 10);
			break;

		default:
			goto usage;
		}
	}

	if (nt < 1) {
		nt = 1;
	}

	Arg *arg = NULL;
	arg = new Arg[nt];
	unsigned c;

	A** ppA = new A*[size];
	for (unsigned i = 0; i < size; i++) {
		ppA[i] = new A(i);
	}

	for (unsigned i = 0; i < nt; i++) {
		arg[i].n = i;
		arg[i].i.Init(ppA, size);
		arg[i].c = &c;
		arg[i].cnt = 0;
		arg[i].stop = false;

		if (i == 0) continue;

		sem_init(&arg[i].s, 0, 0);
		pthread_create(&arg[i].t, NULL, f, &arg[i]);
	}

	for (unsigned k = 0; k < 10; k++) {
		arg[0].i.ResetAccessData();
		c = nt - 1;

		for (unsigned i = 1; i < nt; i++) {
			sem_post(&arg[i].s);
		}

		f2(arg[0].i, arg[0].n, arg[0].cnt);

		if (nt > 1) {
			pthread_mutex_lock(&mutex);
			if (c) {
				silent_cout("main: count " << c << std::endl);
				pthread_cond_wait(&cond, &mutex);
				silent_cout("main: wait is over; count" << c << std::endl);
			}
			pthread_mutex_unlock(&mutex);
		}

		c = arg[0].cnt;
		fprintf(stderr, "cnt = {%u", arg[0].cnt);
		for (unsigned i = 1; i < nt; i++) {
			fprintf(stderr, ", %u", arg[i].cnt);
			c += arg[i].cnt;
		}
		fprintf(stderr, "} = %u\n", c);

		unsigned who[nt];
		memset(who, 0, sizeof(who));
		for (unsigned i = 0; i < size; i++) {
			who[ppA[i]->Get()]++;
		}
	
		c = who[0];
		fprintf(stderr, "who = {%u", who[0]);
		for (unsigned i = 1; i < nt; i++) {
			fprintf(stderr, ", %u", who[i]);
			c += who[i];
		}
		fprintf(stderr, "} = %u\n", c);
	}

	for (unsigned i = 1; i < nt; i++) {
		arg[i].stop = true;
		sem_post(&arg[i].s);
		pthread_join(arg[i].t, NULL);
	}

	delete[] arg;

	return 0;
}

#else /* ! USE_MULTITHREAD */

int
main(void)
{
	std::cerr << "need --enable-multithread" << std::endl;
	exit(EXIT_FAILURE);
}

#endif /* ! USE_MULTITHREAD */

