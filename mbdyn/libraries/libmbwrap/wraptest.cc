/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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

extern "C" {
#include <time.h>
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif /* HAVE_SYS_TIMES_H */
}

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
#include <parsuperluwrap.h>
#include <superluwrap.h>
#include <lapackwrap.h>
#include <taucswrap.h>
#include <naivewrap.h>
#include <parnaivewrap.h>
#include <wsmpwrap.h>

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
#if defined(USE_WSMP)
		"wsmp",
#endif
		"naive",
		NULL
};

std::vector<doublereal> x_values;
std::vector<integer> row_values;
std::vector<integer> col_values;
std::vector<integer> acol_values;
SparseMatrixHandler *spM = 0;

void SetupSystem(
	const bool random,
	char *random_args,
	const bool singular,
	const char *const matrixfilename,
	const char *const filename,
	MatrixHandler *const pM, 
	VectorHandler *const pV
) {
	VariableSubMatrixHandler SBMH(10, 10);
	FullSubMatrixHandler& WM = SBMH.SetFull();
	
	std::ifstream ifile;
	std::ofstream ofile;
	int size = 3;
	
	if (filename == 0) {
		if (random) {
			int size = (*pM).iGetNumRows();
			if (spM == 0) {
				spM = new SpMapMatrixHandler(size);

				int halfband = 0;
				if (random_args) {
					char	*next;

					halfband = (int)strtol(random_args,
							&next, 10);
					switch (next[0]) {
					case '\0':
						random_args = 0;
						break;

					case ':':
						random_args = &next[1];
						break;

					default:
						std::cerr << "unable to parse "
							"<halfband> "
							"from -r args"
							<< std::endl;
						exit(EXIT_FAILURE);
					}

				} else {
					std::cout << "Halfband?" << std::endl;
					std::cin >> halfband;
				}

				int activcol = 0;
				if (random_args) {
					char	*next;

					activcol = (int)strtol(random_args,
							&next, 10);
					switch (next[0]) {
					case '\0':
						random_args = 0;
						break;

					case ':':
						random_args = &next[1];
						break;

					default:
						std::cerr << "unable to parse "
							"<activcol> "
							"from -r args"
							<< std::endl;
						exit(EXIT_FAILURE);
					}

				} else {
					std::cout << "Activcol?" << std::endl;
					std::cin >> activcol;
				}
			
				double sprfct = 0.9;
				if (random_args) {
					char	*next;

					sprfct = (int)strtod(random_args,
							&next);
					switch (next[0]) {
					case '\0':
						random_args = 0;
						break;

					case ':':
						random_args = &next[1];
						break;

					default:
						std::cerr << "unable to parse "
							"<sprfct> "
							"from -r args"
							<< std::endl;
						exit(EXIT_FAILURE);
					}

				} else {
					std::cout << "Sprfct (hint: 0.9)?"
						<< std::endl;
					std::cin >> sprfct;
				}
			
				for (int i = 0; i < size; i++) {
					for (int k = (i - halfband) < 0 ? 0 : i - halfband; k < ((i + halfband) > size ? size : i + halfband); k++) {
						if (((doublereal)rand())/RAND_MAX > sprfct) {
							(*spM)(i+1, k+1) = 2.0*(((doublereal)rand())/RAND_MAX - 0.5);
						}
					}
				}
				for (int i = size - activcol; i < size; i++) {
					for (int k = 0; k < size; k++) {
						if (((doublereal)rand())/RAND_MAX > sprfct) {
							(*spM)(k+1, i+1) = (*spM)(i+1, k+1) = 2.0*(((doublereal)rand())/RAND_MAX - 0.5);
						}
					}
				}
				for (int i = 0; i < size; i++) {
					(*spM)(i+1, i+1) = 1;
				}
				(*spM).MakeIndexForm(x_values, row_values, col_values, acol_values, 1);
				if (matrixfilename != 0) {
					ofile.open(matrixfilename);
					ofile << size << std::endl;
					int n = x_values.size();
					for (int i=0; i<n; i++) {
						ofile << row_values[i] << " " <<  
							col_values[i] << " " << 
							x_values[i] << std::endl;;
					}
					ofile.close();
				}
			}
			int n = x_values.size();
			for (int i=0; i<n; i++) {
				(*pM)(row_values[i], col_values[i]) = x_values[i];
			}
			for (int i = 1; i <= size; i++) {
				pV->PutCoef(i, pM->dGetCoef(i,size));
			}
		} else {

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
		}
	} else {
		ifile.open(filename);
		ifile >> size;

		integer count = 0;
		integer row, col;
		doublereal x;
		while (ifile >> row >> col >> x) {
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
		
		ifile.close();
		
		for (integer i = 1; i <= size; i++) {
			pV->PutCoef(i, pM->dGetCoef(i,size));
		}
		std::cout << "\nThe matrix has "
			<< pM->iGetNumRows() << " rows, "
			<< pM->iGetNumCols() << " cols "
			<< "and " << count << " nonzeros\n" << std::endl;
	}
}

static inline unsigned long long rd_CPU_ts(void)
{
	unsigned long long time;
	__asm__ __volatile__( "rdtsc" : "=A" (time));
	return time;
}

static void
usage(int err)
{
	std::cerr << "usage: wraptest "
		"[-c] "
		"[-d] "
		"[-f <filename>] "
		"[-m <solver>] "
		<< std::endl << "\t\t"
		"[-o] "
		"[-O <option>[=<value>]] "
		"[-p <pivot>] "
		<< std::endl << "\t\t"
		"[-r[<size>[:<halfband>[:<activcol>[:<sprfct>]]]]] "
		<< std::endl << "\t\t"
		"[-s] "
		"[-t <nthreads>] "
		"[-w <filename>] "
		<< std::endl;
	std::cerr << "  -c :  if possible, use compressed column matrix format" << std::endl;
	std::cerr << "\tfor the naive solver: use colamd" << std::endl;
	std::cerr << "  -d :  if possible, use dir matrix format" << std::endl;
	std::cerr << "  -f <filename> : load the matrix from <filename>" << std::endl;
	std::cerr << "\tIf the matrix is loaded from file the solution should be [0 0 .... 1]" << std::endl;
	std::cerr << "\tThe file format is: size row col x row col x etc..." << std::endl;
	std::cerr << "\tThe file format is: size row col x row col x etc..." << std::endl;
	std::cerr << "  -m <solver> : {" << solvers[0];
	for (int i = 1; solvers[i]; i++) {
		std::cerr << "|" << solvers[i];
	}
	std::cerr << "}" << std::endl;
	std::cerr << "  -o :  output of the solution" << std::endl;
	std::cerr << "  -O <option[=<value>]>" << std::endl
		<< "\tblocksize=<blocksize> (umfpack only)" << std::endl;
	std::cerr << "  -p <pivot> : if meaningful, use <pivot> thresold" << std::endl;
	std::cerr << "  -r[<size>[:<halfband>[:<activcol>[:<sprfct>]]]] :" << std::endl;
	std::cerr << "\tgenerate a random matrix with <size>, <halfband>, <activcol>" << std::endl;
	std::cerr << "\tand <sprfct> (prompts for values not provided)" << std::endl;
	std::cerr << "  -s :  (singular) with the 3x3 matrix, do not set the element (3,3)" << std::endl;
	std::cerr << "  -t :  with multi-threaded solutors, use <nthreads> threads" << std::endl;
	std::cerr << "  -w :  write the random matrix to <filename>" << std::endl;
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
#elif defined(USE_WSMP)
		"wsmp"
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
	char *matrixfilename = 0;
	std::ifstream file;
	bool random(false);
	char *random_args = 0;
	bool cc(false);
	bool dir(false);
	unsigned nt = 1;
	unsigned block_size = 0;
	double dpivot = -1.;
	bool singular(false);
	bool output_solution(false);
	int size = 3;
	long long tf;
	
	while (1) {
		int opt = getopt(argc, argv, "cdf:m:oO:p:r::st:w:");

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

		case 'p':
			dpivot = atof(optarg);
			break;

		case 'O':
		{
			if (strncasecmp(optarg, "blocksize=", sizeof("blocksize=") - 1) == 0) {
				char	*next;

				if (strcasecmp(solver, "umfpack") != 0) {
					std::cerr << "blocksize only meaningful for umfpack solver" << std::endl;
				}

				optarg += sizeof("blocksize=") - 1;
				block_size = (int)strtol(optarg, &next, 10);
				if (next[0] != '\0') {
					std::cerr << "unable to parse blocksize value" << std::endl;
					exit(EXIT_FAILURE);
				}
				if (block_size < 1) {
					block_size = 0;
				}

			} else {
				std::cerr << "unrecognized option \"" << optarg << "\"" << std::endl;
				exit(EXIT_FAILURE);
			}

			break;
		}

		case 'f':
			filename = optarg;
			break;

		case 'o':
			output_solution = true;
			break;

		case 'r':
			random = true;
			random_args = optarg;
			break;
			
		case 'w':
			matrixfilename = optarg;
			break;


		default:
			usage(EXIT_FAILURE);
		}
	}

	if (filename != 0) {
		file.open(filename);
		file >> size;
		file.close();

	} else if (random) {
		if (random_args) {
			char	*next;

			size = (int)strtol(random_args, &next, 10);
			switch (next[0]) {
			case '\0':
				random_args = 0;
				break;

			case ':':
				random_args = &next[1];
				break;

			default:
				std::cerr << "unable to parse <size> "
					"from -r args" << std::endl;
				exit(EXIT_FAILURE);
			}

		} else {
			std::cout << "Matrix size?" << std::endl;
			std::cin >> size;
		}
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
		if (nt > 1) {
#ifdef USE_SUPERLU_MT
			std::cerr << "Multi-threaded SuperLU solver";
			if (dir) {
				std::cerr << " with dir matrix";
				typedef ParSuperLUSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
				SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(nt, size));
			} else if (cc) {
				std::cerr << " with cc matrix";
				typedef ParSuperLUSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
				SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(nt, size));
			} else {
				SAFENEWWITHCONSTRUCTOR(pSM, ParSuperLUSparseSolutionManager,
						ParSuperLUSparseSolutionManager(nt, size));
			}
			std::cerr << " using " << nt << " threads" << std::endl;
#else /* !USE_SUPERLU_MT */
			silent_cerr("multithread SuperLU solver support not compiled; "
				<< std::endl);
			throw ErrGeneric();
		} else {
			std::cerr << "SuperLU solver";
			if (dir) {
				std::cerr << " with dir matrix";
				typedef SuperLUSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
				SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));
			} else if (cc) {
				std::cerr << " with cc matrix";
				typedef SuperLUSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
				SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));
			} else {
				SAFENEWWITHCONSTRUCTOR(pSM, SuperLUSparseSolutionManager,
					SuperLUSparseSolutionManager(size));
			}
#endif /* !USE_SUPERLU_MT */
		}
		std::cerr << std::endl;
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
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size, dpivot, block_size));

		} else if (cc) {
			std::cerr << " with cc matrix";
			typedef UmfpackSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size, dpivot, block_size));

		} else {
			SAFENEWWITHCONSTRUCTOR(pSM,
					UmfpackSparseSolutionManager,
					UmfpackSparseSolutionManager(size, dpivot, block_size));
		}
		std::cerr << std::endl;
#else /* !USE_UMFPACK */
		std::cerr << "need --with-umfpack to use Umfpack library" 
			<< std::endl;
		usage(EXIT_FAILURE);
#endif /* !USE_UMFPACK */

	} else if (strcasecmp(solver, "wsmp") == 0) {
#ifdef USE_WSMP
		std::cerr << "Wsmp solver";
		if (dir) {
			std::cerr << " with dir matrix";
			typedef WsmpSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size, dpivot, block_size, nt));

		} else if (cc) {
			std::cerr << " with cc matrix";
			typedef WsmpSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size, dpivot, block_size, nt));

		} else {
			SAFENEWWITHCONSTRUCTOR(pSM,
					WsmpSparseSolutionManager,
					WsmpSparseSolutionManager(size, dpivot, block_size, nt));
		}
		std::cerr << " using " << nt << " threads " << std::endl;
		std::cerr << std::endl;
#else /* !USE_WSMP */
		std::cerr << "need --with-wsmp to use Wsmp library" 
			<< std::endl;
		usage(EXIT_FAILURE);
#endif /* !USE_WSMP */

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

	SetupSystem(random, random_args, singular,
			matrixfilename, filename, pM, pV);
	
	clock_t ct = 0;
	struct tms tmsbuf;
	try {
		start = clock();
		tf = rd_CPU_ts();
		times(&tmsbuf);
		ct = tmsbuf.tms_utime + tmsbuf.tms_cutime
			+ tmsbuf.tms_stime + tmsbuf.tms_cstime;
		pSM->Solve();
		tf = rd_CPU_ts() - tf;
		end = clock();
		times(&tmsbuf);
		ct = tmsbuf.tms_utime + tmsbuf.tms_cutime
			+ tmsbuf.tms_stime + tmsbuf.tms_cstime - ct;
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
	std::cout << "Clock tics to solve: " << tf << std::endl;
	std::cout << "Clock tics to solve: " << ct << std::endl;
	std::cout << "Time to solve: " << cpu_time_used << std::endl;
	

	std::cout << "\nSecond solve:\n";
	
	pSM->MatrReset();
	pM = pSM->pMatHdl();

	pM->Reset();

	SetupSystem(random, random_args, singular, 0, filename, pM, pV);
	
	try {
		start = clock();
		tf = rd_CPU_ts();
		times(&tmsbuf);
		ct = tmsbuf.tms_utime + tmsbuf.tms_cutime
			+ tmsbuf.tms_stime + tmsbuf.tms_cstime;
		pSM->Solve();
		tf = rd_CPU_ts() - tf;
		end = clock();
		times(&tmsbuf);
		ct = tmsbuf.tms_utime + tmsbuf.tms_cutime
			+ tmsbuf.tms_stime + tmsbuf.tms_cstime - ct;
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
	std::cout << "Clock tics to solve: " << tf << std::endl;
	std::cout << "Clock tics to solve: " << ct << std::endl;
	std::cout << "Time to solve: " << cpu_time_used << std::endl;
	
	return 0;
}

