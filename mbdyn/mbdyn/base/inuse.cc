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

void *
f(void *p)
{
	Arg *arg = (Arg *)p;

	sem_wait(&arg->s);

	A* pA = NULL;
	if (arg->i.bGetFirst(pA)) {
		do {
			pA->Set(arg->n);
		} while (arg->i.bGetNext(pA));
	}

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
	unsigned c = nt;

	A** ppA = new A*[size];
	for (unsigned i = 0; i < size; i++) {
		ppA[i] = new A;
		//ppA[i]->SetInUse(false);
	}

	for (unsigned i = 0; i < nt; i++) {
		sem_init(&arg[i].s, 0, 0);
		arg[i].n = i;
		arg[i].i.Init(ppA, size);
		arg[i].c = &c;

		pthread_create(&arg[i].t, NULL, f, &arg[i]);
	}

	pthread_mutex_lock(&mutex);

	for (unsigned i = 0; i < nt; i++) {
		sem_post(&arg[i].s);
	}

	pthread_cond_wait(&cond, &mutex);
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

