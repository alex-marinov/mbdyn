/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
#include <mbconfig.h>
#endif /* HAVE_CONFIG */

#include <stdlib.h>
#include <ac/iostream>

#include <solman.h>
#include <y12wrap.h>
#include <harwrap.h>
#include <mschwrap.h>
#include <umfpackwrap.h>

static void
usage(void)
{
	std::cerr << "usage: wraptest "
		"[y12|harwell|meschach|{umfpack [cc]} [singular]]" 
		<< std::endl;
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
		SAFENEWWITHCONSTRUCTOR(pSM, Y12SparseLUSolutionManager,
				Y12SparseLUSolutionManager(size));
#else /* !USE_Y12 */
		std::cerr << "need --with-y12 to use y12m library" 
			<< std::endl;
		usage();
#endif /* !USE_Y12 */

	} else if (strcasecmp(solver, "harwell") == 0) {
#ifdef USE_HARWELL
		SAFENEWWITHCONSTRUCTOR(pSM, HarwellSparseLUSolutionManager,
				HarwellSparseLUSolutionManager(size));
#else /* !USE_HARWELL */
		std::cerr << "need --with-harwell to use HSL library" 
			<< std::endl;
		usage();
#endif /* !USE_HARWELL */

	} else if (strcasecmp(solver, "meschach") == 0) {
#ifdef USE_MESCHACH
		SAFENEWWITHCONSTRUCTOR(pSM, MeschachSparseLUSolutionManager,
				MeschachSparseLUSolutionManager(size));
#else /* !USE_MESCHACH */
		std::cerr << "need --with-meschach to use Meschach library" 
			<< std::endl;
		usage();
#endif /* !USE_MESCHACH */
	} else if (strcasecmp(solver, "umfpack") == 0
			|| strcasecmp(solver, "umfpack3") == 0) {
#ifdef USE_UMFPACK
		if (argc > 2 && strcasecmp(argv[2], "cc") == 0) {
			SAFENEWWITHCONSTRUCTOR(pSM,
					UmfpackSparseCCLUSolutionManager,
					UmfpackSparseCCLUSolutionManager(size));
			argc--;
			argv++;
		} else {
			SAFENEWWITHCONSTRUCTOR(pSM,
					UmfpackSparseLUSolutionManag,
					UmfpackSparseLUSolutionManager(size));
		}
#else /* !USE_UMFPACK */
		std::cerr << "need --with-umfpack to use Umfpack library" 
			<< std::endl;
		usage();
#endif /* !USE_UMFPACK */

	} else {
		std::cerr << "unknown solver '" << solver << "'" << std::endl;
		usage();
	}
	argc--;
	argv++;

	pSM->MatrInit(0.);
	
	MatrixHandler *pM = pSM->pMatHdl();
	VectorHandler *pV = pSM->pResHdl();
	VectorHandler *px = pSM->pSolHdl();

	pM->PutCoef(1, 1, 1.);
	pM->PutCoef(2, 2, 2.);
	if (argc > 1 && strcasecmp(argv[1], "singular") == 0) {
		pM->PutCoef(3, 3, 0.);
	} else {
		pM->PutCoef(3, 3, 3.);
	}
	
	pV->PutCoef(1, 1.);
	pV->PutCoef(2, 1.);
	pV->PutCoef(3, 1.);
	
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
		std::cout << "\tsol[" << i << "] = " << px->dGetCoef(i) 
			<< std::endl;
	}
	
	return 0;
}

