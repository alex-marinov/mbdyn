/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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
#include <unistd.h>

#include <time.h>

#include <ac/iostream>
#include <ac/fstream>

#include <solman.h>
#include <spmapmh.h>
#include <ccmh.h>
#include <dirccmh.h>
#include <y12wrap.h>
#include <harwrap.h>
#include <mschwrap.h>
#include <umfpackwrap.h>
#include <superluwrap.h>
#include <lapackwrap.h>
#include <taucswrap.h>
#include <naivewrap.h>
#include <parnaivewrap.h>

char *solvers[] = {
#if defined(USE_Y12)
		"y12",
#endif
#if defined(USE_UMFPACK)
		"umfpack",
#endif
#if defined(USE_SUPERLU)
		"superlu",
#endif
#if defined(USE_HARWELL)
		"harwell",
#endif
#if defined(USE_MESCHACH)
		"meschach",
#endif
#if defined(USE_LAPACK)
		"lapack",
#endif
#if defined(USE_TAUCS)
		"taucs",
#endif
		"naive",
		NULL
};

void SetupSystem(
	const bool singular,
	const char *const filename,
	MatrixHandler *const pM, 
	VectorHandler *const pV
) {
	VariableSubMatrixHandler SBMH(10, 10);
	FullSubMatrixHandler& WM = SBMH.SetFull();
	
	std::ifstream file;
	int size(3);
	
	if (filename == 0) {

		WM.ResizeReset(3, 3);
		WM.PutRowIndex(1,1);
		WM.PutRowIndex(2,2);
		WM.PutRowIndex(3,3);
		WM.PutColIndex(1,1);
		WM.PutColIndex(2,2);
		WM.PutColIndex(3,3);
		WM.PutCoef(1, 1, 1.);
		WM.PutCoef(2, 2, 2.);
		WM.PutCoef(2, 3, -1.);
		WM.PutCoef(3, 2, 11.);
		WM.PutCoef(3, 1, 10.);
		if (singular) {
			WM.PutCoef(3, 3, 0.);

		} else {
			WM.PutCoef(3, 3, 3.);
		}
	
		*pM += SBMH;
	
		pV->PutCoef(1, 1.);
		pV->PutCoef(2, 1.);
		pV->PutCoef(3, 1.);
	
		std::cout << *pM << "\n";
	} else {
		file.open(filename);
		file >> size;

		integer count = 0;
		integer row, col;
		doublereal x;
		while (file >> row >> col >> x) {
			if (row > size || col > size) {
				std::cerr << "Fatal read error of file" << filename << std::endl;
				std::cerr << "size: " << size << std::endl;
				std::cerr << "row:  " << row << std::endl;
				std::cerr << "col:  " << col << std::endl;
				std::cerr << "x:    " << x << std::endl;
			}
			(*pM)(row, col) = x;
			count++;
		}
		
		file.close();
		
		for (integer i = 1; i <= size; i++) {
			pV->PutCoef(i, (*pM)(i,size));
		}
		std::cout << "\nThe matrix has "
			<< pM->iGetNumRows() << " rows, "
			<< pM->iGetNumCols() << " cols "
			<< "and " << count << " nonzeros\n" << std::endl;
	}
}

// static inline unsigned long long rd_CPU_ts(void)
// {
// 	unsigned long long time;
// 	__asm__ __volatile__( "rdtsc" : "=A" (time));
// 	return time;
// }

static void
usage(int err)
{
	std::cerr << "usage: wraptest [-c] [-d] [-m <solver>] [-s] [-t <nthreads>] [-f <filename>] [-o]" << std::endl
		<< "\t<solver>={" << solvers[0];
	std::cerr << "If the matrix is loaded from file the solution should be [0 0 .... 1]" << std::endl;
	std::cerr << "The file format is: size row col x row col x etc..." << std::endl;
	for (int i = 1; solvers[i]; i++) {
		std::cerr << "|" << solvers[i];
	}
	std::cerr << "}" << std::endl;
	exit(err);
}

int
main(int argc, char *argv[])
{
	SolutionManager *pSM = NULL;
	char *solver =
#if defined(USE_UMFPACK)
		"umfpack"
#elif defined(USE_Y12)
		"y12"
#elif defined(USE_SUPERLU)
		"superlu"
#elif defined(USE_HARWELL)
		"harwell"
#elif defined(USE_MESCHACH)
		"meschach"
#elif defined(USE_LAPACK)
		"lapack"
#elif defined(USE_TAUCS)
		"taucs"
#else
		"naive"
#if 0
		"no solver!!!"
#endif
#endif /* NO SOLVER !!! */
		;
	clock_t start, end;
	double cpu_time_used;
	
	char *filename = 0;
	std::ifstream file;
	bool cc(false);
	bool dir(false);
	unsigned nt = 1;
	bool singular(false);
	bool output_solution(false);
	int size(3);
	long long tf;

	while (1) {
		int opt = getopt(argc, argv, "cdm:st:f:o");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'c':
			cc = true;
			break;

		case 'd':
			dir = true;
			break;

		case 'm':
			solver = optarg;
			break;

		case 's':
			singular = true;
			break;

		case 't':
			nt = atoi(optarg);
			if (nt < 1) {
				nt = 1;
			}
			break;

		case 'f':
			filename = optarg;
			break;

		case 'o':
			output_solution = true;
			break;

		default:
			usage(EXIT_FAILURE);
		}
	}

	if (filename != 0) {
		file.open(filename);
		file >> size;
		file.close();
	}
	
	std::cerr << std::endl;

	if (strcasecmp(solver, "taucs") == 0) {
#ifdef USE_TAUCS
		std::cerr << "Taucs solver"
		if (dir) {
			std::cerr << " with dir matrix"
			typedef TaucsSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else if (cc) {
			std::cerr << " with cc matrix"
			typedef TaucsSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else {
			SAFENEWWITHCONSTRUCTOR(pSM, TaucsSparseSolutionManager,
					TaucsSparseSolutionManager(size));
		}
		std::cerr << std::endl;
#else /* !USE_TAUCS */
		std::cerr << "need --with-taucs to use Taucs library sparse solver" 
			<< std::endl;
		usage(EXIT_FAILURE);
#endif /* !USE_LAPACK */

	} else if (strcasecmp(solver, "lapack") == 0) {
#ifdef USE_LAPACK
		std::cerr << "Lapack solver" << std::endl;
		SAFENEWWITHCONSTRUCTOR(pSM, LapackSolutionManager,
				LapackSolutionManager(size));
#else /* !USE_LAPACK */
		std::cerr << "need --with-lapack to use Lapack library dense solver" 
			<< std::endl;
		usage(EXIT_FAILURE);
#endif /* !USE_LAPACK */

	} else if (strcasecmp(solver, "superlu") == 0) {
#ifdef USE_SUPERLU
		std::cerr << "SuperLU solver";
		if (dir) {
			std::cerr << " with dir matrix";
			typedef SuperLUSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(nt, size));

		} else if (cc) {
			std::cerr << " with cc matrix";
			typedef SuperLUSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(nt, size));

		} else {
			SAFENEWWITHCONSTRUCTOR(pSM, SuperLUSparseSolutionManager,
					SuperLUSparseSolutionManager(nt, size));
		}
		std::cerr << " using " << nt << " threads" << std::endl;
#else /* !USE_SUPERLU */
		std::cerr << "need --with-superlu to use SuperLU library" 
			<< std::endl;
		usage(EXIT_FAILURE);
#endif /* !USE_SUPERLU */

	} else if (strcasecmp(solver, "y12") == 0) {
#ifdef USE_Y12
		std::cerr << "y12 solver";
		if (dir) {
			std::cerr << " with dir matrix";
			typedef Y12SparseCCSolutionManager<DirCColMatrixHandler<1> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else if (cc) {
			std::cerr << " with cc matrix";
			typedef Y12SparseCCSolutionManager<CColMatrixHandler<1> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else {
			SAFENEWWITHCONSTRUCTOR(pSM, Y12SparseSolutionManager,
					Y12SparseSolutionManager(size));
		}
		std::cerr << std::endl;
#else /* !USE_Y12 */
		std::cerr << "need --with-y12 to use y12m library" 
			<< std::endl;
		usage(EXIT_FAILURE);
#endif /* !USE_Y12 */

	} else if (strcasecmp(solver, "harwell") == 0) {
#ifdef USE_HARWELL
		std::cerr << "Harwell solver" << std::endl;
		SAFENEWWITHCONSTRUCTOR(pSM, HarwellSparseSolutionManager,
				HarwellSparseSolutionManager(size));
#else /* !USE_HARWELL */
		std::cerr << "need --with-harwell to use HSL library" 
			<< std::endl;
		usage(EXIT_FAILURE);
#endif /* !USE_HARWELL */

	} else if (strcasecmp(solver, "meschach") == 0) {
#ifdef USE_MESCHACH
		std::cerr << "Meschach solver" << std::endl;
		SAFENEWWITHCONSTRUCTOR(pSM, MeschachSparseSolutionManager,
				MeschachSparseSolutionManager(size));
#else /* !USE_MESCHACH */
		std::cerr << "need --with-meschach to use Meschach library" 
			<< std::endl;
		usage(EXIT_FAILURE);
#endif /* !USE_MESCHACH */
	} else if (strcasecmp(solver, "umfpack") == 0
			|| strcasecmp(solver, "umfpack3") == 0) {
#ifdef USE_UMFPACK
		std::cerr << "Umfpack solver";
		if (dir) {
			std::cerr << " with dir matrix";
			typedef UmfpackSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else if (cc) {
			std::cerr << " with cc matrix";
			typedef UmfpackSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else {
			SAFENEWWITHCONSTRUCTOR(pSM,
					UmfpackSparseSolutionManager,
					UmfpackSparseSolutionManager(size));
		}
		std::cerr << std::endl;
#else /* !USE_UMFPACK */
		std::cerr << "need --with-umfpack to use Umfpack library" 
			<< std::endl;
		usage(EXIT_FAILURE);
#endif /* !USE_UMFPACK */

	} else if (strcasecmp(solver, "naive") == 0) {
		std::cerr << "Naive solver";
		if (cc) {
			std::cerr << " with Colamd ordering";
			if (nt > 1) {
#ifdef USE_NAIVE_MULTITHREAD
				SAFENEWWITHCONSTRUCTOR(pSM,
					ParNaiveSparsePermSolutionManager,
					ParNaiveSparsePermSolutionManager(nt, size, 1.E-8));
#else
				silent_cerr("multithread naive solver support not compiled; "
					"you can configure --enable-multithread-naive "
					"on a linux ix86 to get it"
					<< std::endl);
				throw ErrGeneric();
#endif /* USE_NAIVE_MULTITHREAD */
			} else {
				SAFENEWWITHCONSTRUCTOR(pSM,
					NaiveSparsePermSolutionManager,
					NaiveSparsePermSolutionManager(size, 1.E-8));
			}
		} else {
			if (nt > 1) {
#ifdef USE_NAIVE_MULTITHREAD
				SAFENEWWITHCONSTRUCTOR(pSM,
					ParNaiveSparseSolutionManager,
					ParNaiveSparseSolutionManager(nt, size, 1.E-8));
#else
				silent_cerr("multithread naive solver support not compiled; "
					"you can configure --enable-multithread-naive "
					"on a linux ix86 to get it"
					<< std::endl);
				throw ErrGeneric();
#endif /* USE_NAIVE_MULTITHREAD */
			} else {
				SAFENEWWITHCONSTRUCTOR(pSM,
					NaiveSparseSolutionManager,
					NaiveSparseSolutionManager(size, 1.E-8));
			}
		}
		std::cerr << " using " << nt << " threads " << std::endl;
	} else {
		std::cerr << "unknown solver '" << solver << "'" << std::endl;
		usage(EXIT_FAILURE);
	}

	pSM->MatrReset();
	
	MatrixHandler *pM = pSM->pMatHdl();
	VectorHandler *pV = pSM->pResHdl();
	VectorHandler *px = pSM->pSolHdl();

	pM->Reset();

	SetupSystem(singular, filename, pM, pV);
	
	try {
		start = clock();
		//tf = rd_CPU_ts();
		pSM->Solve();
		//tf = rd_CPU_ts() - tf;
		end = clock();
	} catch (...) {
		exit(EXIT_FAILURE);
	}
	if (output_solution) {
		for (int i = 1; i <= size; i++) {
			std::cout << "\tsol[" << i << "] = " << px->dGetCoef(i) 
				<< std::endl;
		}
	}
	//cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cpu_time_used = ((double) (end - start));
	cpu_time_used = cpu_time_used / CLOCKS_PER_SEC;
	std::cout << "Clock tics to solve: " << end - start << std::endl;
	std::cout << "Time to solve: " << cpu_time_used << std::endl;
	

	std::cout << "\nSecond solve:\n";
	
	pSM->MatrReset();
	pM = pSM->pMatHdl();

	pM->Reset();

	SetupSystem(singular, filename, pM, pV);
	
	try {
		start = clock();
		//tf = rd_CPU_ts();
		pSM->Solve();
		//tf = rd_CPU_ts() - tf;
		end = clock();
	} catch (...) {
		exit(EXIT_FAILURE);
	}
	
	if (output_solution) {
		for (int i = 1; i <= size; i++) {
			std::cout << "\tsol[" << i << "] = " << px->dGetCoef(i) 
				<< std::endl;
		}
	}
	//cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cpu_time_used = ((double) (end - start));
	cpu_time_used = cpu_time_used / CLOCKS_PER_SEC;
	std::cout << "Clock tics to solve: " << end - start << std::endl;
	std::cout << "Time to solve: " << cpu_time_used << std::endl;
	
	return 0;
}

