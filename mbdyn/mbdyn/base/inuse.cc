/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <ac/iostream>

#ifdef USE_MULTITHREAD

#include <unistd.h>
#include <pthread.h>
#include <semaphore.h>

#include "mbconfig.h"
#include "veciter.h"

unsigned dst = 100;
unsigned rst = 10;

class A : public InUse {
private:
	unsigned who;
public:
	A(void) {};
	virtual ~A(void) {};

	void Set(unsigned w) { who = w; usleep(dst + rand()%rst); };
	unsigned Get(void) const { return who; };
};

struct Arg {
	unsigned n;
	pthread_t t;
	sem_t s;
	VecIter<A *> i;
	unsigned *c;
};

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;

static void
f2(VecIter<A *>* pi, unsigned n)
{
	A* pA = NULL;

	if (pi->bGetFirst(pA)) {
		do {
			pA->Set(n);
		} while (pi->bGetNext(pA));
	}
}

void *
f(void *p)
{
	Arg *arg = (Arg *)p;

	sem_wait(&arg->s);

	f2(&arg->i, arg->n);

	pthread_mutex_lock(&mutex);
	--*(arg->c);
	if (*arg->c == 0) {
		pthread_cond_signal(&cond);
	}
	pthread_mutex_unlock(&mutex);

	return NULL;
}

int
main(int argc, char* argv[])
{
	unsigned size = 1000;
	unsigned nt = 1;

	if (argc > 1) {
		size = atoi(argv[1]);
		if (size <= 0) {
			return -1;
		}
		
		if (argc > 2) {
			nt = atoi(argv[2]);
			if (nt <= 0) {
				return -1;
			}

			if (argc > 3) {
				dst = atoi(argv[3]);
				if (dst <= 0) {
					return -1;
				}

				if (argc > 4) {
					rst = atoi(argv[4]);
					if (rst <= 0) {
						return -1;
					}
				}
			}
		}
	}

	Arg arg[nt];
	unsigned c = nt - 1;

	A** ppA = new A*[size];
	for (unsigned i = 0; i < size; i++) {
		ppA[i] = new A;
		ppA[i]->SetInUse(false);
	}

	for (unsigned i = 0; i < nt; i++) {
		arg[i].n = i;
		arg[i].i.Init(ppA, size);
		arg[i].c = &c;

		if (i == 0) continue;

		sem_init(&arg[i].s, 0, 0);
		pthread_create(&arg[i].t, NULL, f, &arg[i]);
	}

	pthread_mutex_lock(&mutex);

	for (unsigned i = 1; i < nt; i++) {
		sem_post(&arg[i].s);
	}

	f2(&arg[0].i, arg[0].n);

	if (nt > 1) {
		pthread_cond_wait(&cond, &mutex);
	}
	pthread_mutex_unlock(&mutex);

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

