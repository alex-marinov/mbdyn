#ifdef HAVE_CONFIG_H
#include <mbconfig.h>
#endif

#include <stdlib.h>
#include <iostream.h>
#include <string.h>

#include <solman.h>
#include <y12wrap.h>
#include <harwrap.h>
#include <mschwrap.h>

static void
usage(void)
{
	cerr << "usage: t [y12|harwell|meschach [singular]]" << endl;
	exit(EXIT_FAILURE);
}

int
main(int argc, char *argv[])
{
	SolutionManager *pSM = NULL;
	char *solver = "y12";
	const int size = 3;

	if (argc > 1) {
		solver = argv[1];
	}

	if (strcasecmp(solver, "y12") == 0) {
#ifdef USE_Y12
		pSM = new Y12SparseLUSolutionManager(size);
#else 
		cerr << "need --with-y12 to use y12m library" << endl;
		usage();
#endif

	} else if (strcasecmp(solver, "harwell") == 0) {
#ifdef USE_HARWELL
		pSM = new HarwellSparseLUSolutionManager(size);
#else
		cerr << "need --with-harwell to use HSL library" << endl;
		usage();
#endif

	} else if (strcasecmp(solver, "meschach") == 0) {
#ifdef USE_MESCHACH
		pSM = new MeschachSparseLUSolutionManager(size);
#else
		cerr << "need --with-meschach to use Meschach library" << endl;
		usage();
#endif
	} else {
		cerr << "unknown solver '" << solver << "'" << endl;
		usage();
	}

	pSM->MatrInit(0.);
	
	MatrixHandler *pM = pSM->pMatHdl();
	VectorHandler *pV = pSM->pResHdl();

	pM->fPutCoef(1, 1, 1.);
	pM->fPutCoef(2, 2, 2.);
	if (argc > 2 && strcasecmp(argv[2], "singular") == 0) {
		pM->fPutCoef(3, 3, 0.);
	} else {
		pM->fPutCoef(3, 3, 3.);
	}
	
	pV->fPutCoef(1, 1.);
	pV->fPutCoef(2, 1.);
	pV->fPutCoef(3, 1.);
	
#ifdef USE_EXCEPTIONS
	try {
#endif /* USE_EXCEPTIONS */
	pSM->Solve();
#ifdef USE_EXCEPTIONS
	} catch (...) {
		exit(EXIT_FAILURE);
	}
#endif /* USE_EXCEPTIONS */
	
	for (int i = 1; i <= size; i++) {
		cout << "\tsol[" << i << "] = " << pV->dGetCoef(i) << endl;
	}
	
	return 0;
}

