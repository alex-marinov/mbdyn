/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "mbconfig.h"

#include <cstring>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include "ac/getopt.h"

#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "y12wrap.h"
#include "harwrap.h"
#include "mschwrap.h"
#include "umfpackwrap.h"
#include "parsuperluwrap.h"
#include "superluwrap.h"

const char *solvers[] = {
#ifdef USE_UMFPACK
	"umfpack",
#endif /* USE_UMFPACK */
#ifdef USE_Y12
	"y12",
#endif /* USE_Y12 */
#ifdef USE_SUPERLU
	"superlu",
#endif /* USE_SUPERLU */
	NULL
};

static void
usage(void)
{
	std::cerr << "usage: cctest [-d] [-m <solver>] [-t <nthreads>]"
		<< std::endl;

	if (!solvers[0]) {
		std::cerr << "\tno solvers available!!!" << std::endl;

	} else {
		std::cerr << "\t<solver>={" << solvers[0];
		for (unsigned i = 1; solvers[i]; i++) {
			std::cerr << "|" << solvers[i];
		}
		std::cerr << "}" << std::endl;
	}

	exit(EXIT_FAILURE);
}

int
main(int argc, char *argv[])
{
	SolutionManager *pSM = NULL;
	char *solver = NULL;
	bool dir(false);
	unsigned nt = 1;
	const int size(3);

	while (1) {
		int opt = getopt(argc, argv, "dm:t:");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'd':
			dir = true;
			break;

		case 'm':
			solver = optarg;
			break;

		case 't':
			nt = atoi(optarg);
			if (nt < 1) {
				nt = 1;
			}
			break;

		default:
			usage();
		}
	}

	if (solver == NULL) {
		usage();
	}

	if (strcasecmp(solver, "superlu") == 0) {
#ifdef USE_SUPERLU
#ifdef USE_SUPERLU_MT
		if (dir) {
			typedef ParSuperLUSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(nt, size));

		} else {
			typedef ParSuperLUSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(nt, size));
		}
#else /* !USE_SUPERLU_MT */
		if (dir) {
			typedef SuperLUSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else {
			typedef SuperLUSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));
		}
#endif /* USE_SUPERLU_MT */
#else /* !USE_SUPERLU */
		std::cerr << "need --with-superlu to use SuperLU library" 
			<< std::endl;
		usage();
#endif /* !USE_SUPERLU */

	} else if (strcasecmp(solver, "y12") == 0) {
#ifdef USE_Y12
		if (dir) {
			typedef Y12SparseCCSolutionManager<DirCColMatrixHandler<1> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else {
			typedef Y12SparseCCSolutionManager<CColMatrixHandler<1> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));
		}
#else /* !USE_Y12 */
		std::cerr << "need --with-y12 to use y12m library" 
			<< std::endl;
		usage();
#endif /* !USE_Y12 */

	} else if (strcasecmp(solver, "umfpack") == 0
			|| strcasecmp(solver, "umfpack3") == 0) {
#ifdef USE_UMFPACK
		if (dir) {
			typedef UmfpackSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else {
			typedef UmfpackSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));
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

	std::cout << "using " << solver << " solver..." << std::endl;

	MatrixHandler *pM = pSM->pMatHdl();
	VectorHandler *pV = pSM->pResHdl();
	VectorHandler *px = pSM->pSolHdl();

	int count = 0;

	int p = 2;
	int w = 7 + p;
	std::cout.setf(std::ios::scientific);
	std::cout.precision(p);

	double	m11 = 4.,
		m22 = 4.,
		m33 = 2.,
		d = 0.,
		b1 = 0.,
		b2 = 0.,
		b3 = 1.;

retry:;
	pSM->MatrReset();

	try {
		pM = pSM->pMatHdl();
		pM->Reset();

		pM->operator()(1, 1) = m11;
		pM->operator()(2, 2) = m22;
		pM->operator()(3, 3) = m33;
		if (d) {
			pM->operator()(1, 2) = d;
			pM->operator()(2, 1) = d;
			pM->operator()(2, 3) = d;
			pM->operator()(3, 2) = d;
		}
		
		pV->operator()(1) = b1;
		pV->operator()(2) = b2;
		pV->operator()(3) = b3;

	} catch (MatrixHandler::ErrRebuildMatrix) {
		std::cerr << "need to rebuild matrix..." << std::endl;
		pSM->MatrInitialize();
		goto retry;

	} catch (...) {
		std::cerr << "build failure" << std::endl;
		exit(EXIT_FAILURE);
	}

	try {
		count++;
		std::cout << "solution " << count << "..." << std::endl;
		
		pSM->Solve();

	} catch (...) {
		std::cerr << "solution failure" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::cout
		<< "{x1}   [" << std::setw(w) << m11
			<< "," << std::setw(w) << d 
			<< "," << std::setw(w) << 0.
		<< "]^-1 {" << std::setw(w) << b1
		<< "}   {" << std::setw(w) << (*px)(1) << "}" << std::endl
		<< "{x2} = ["<< std::setw(w) << d
			<< "," << std::setw(w) << m22
			<< "," << std::setw(w) << d 
		<< "]    {" << std::setw(w) << b2
		<< "} = {" << std::setw(w) << (*px)(2) << "}" << std::endl	
		<< "{x3}   [" << std::setw(w) << 0.
			<< "," << std::setw(w) << d 
			<< "," << std::setw(w) << m33
		<< "]    {" << std::setw(w) << b3
		<< "}   {" << std::setw(w) << (*px)(3) << "}" << std::endl;

	switch (count) {
	case 1:
		d = -2.;
		break;

	case 2:
		d = 0.;
		break;

	case 3:
		return 0;
	}

	goto retry;
}

