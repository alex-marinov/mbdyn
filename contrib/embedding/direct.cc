/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include "mbconfig.h"

#include "mbdyn.h"
#include "solver.h"

#define BUFFER_STL 0x00U
#define BUFFER_RAW 0x01U
#define BUFFER_MBDYN_OWNS_MEM 0x10U

// use "type, stl"
//#define BUFFER_WHICH BUFFER_STL
//
// use "type, raw, owns memory, no"
#define BUFFER_WHICH BUFFER_RAW
//
// use "type, raw"
//#define BUFFER_WHICH (BUFFER_RAW|BUFFER_MBDYN_OWNS_MEM)

int
main(int argc, char *argv[])
{
	const char *sIn = "direct.mbd";
	const char *sOut = "output";

	Table T(true);
	MathParser MP(T, false);
	std::ifstream streamIn(sIn);
	InputStream In(streamIn);
	MBDynParser HP(MP, In, sIn);
	Solver *pS = new Solver(HP, sIn, sOut, false);

	::fSilent = 0;
	const char *sOutFileName = 0;
	double c = 0.; // extra damping

	while (1) {
		int iCurrOpt = getopt(argc, argv, "c:ho:s");

		if (iCurrOpt == EOF) {
			break;
		}

		switch (iCurrOpt) {
		case int('c'):
			{
				char *endptr;
				c = strtod(optarg, &endptr);
				if (endptr == optarg) {
					std::cerr << "conversion error of damping value '" << optarg << "'" << std::endl;
					exit(1);
				}
			}
			break;
			
		case int('o'):
			sOutFileName = optarg;
			break;

		case int('s'):
			::fSilent++;
			break;

		case int('h'):
		default:
			if (iCurrOpt != int('h')) {
				std::cerr << "unknown option '" << char(iCurrOpt) << "'" << std::endl;
			}
			std::cerr
				<< "usage:" << std::endl
				<< "    -c <damping> # additional damping" << std::endl
				<< "    -o <output file name> # output file name (default: stdout)" << std::endl
				<< "    -s # more silent" << std::endl
				<< std::endl;
			exit(0);
		}
	}

	std::ofstream fout;
	std::ostream *out = &std::cout;
	if (sOutFileName) {
		fout.open(sOutFileName);
		out = &fout;
	}

	if (pS->Prepare()) {
		// connect I/O buffers
		DataManager *pDM = pS->pGetDataManager();

		// define size of I/O buffers, create and currently read them from file

		// I/O stream labels - must be consistent with labels used in the model
		const unsigned uInKinLabel = 97;
		const unsigned uInFrcLabel = 98;
		const unsigned uOutLabel = 0;

		// I/O stream sizes - must be consistent with sizes used in the model
		const int iInKinSize = 1; // prescribed displacement, scalar
		const int iInFrcSize = 1; // prescribed force, scalar
		const int iOutSize = 2; // displacement and velocity


#if BUFFER_WHICH == BUFFER_STL
		std::vector<doublereal>& inKinbuf = pDM->GetBufIn(uInKinLabel);
		std::vector<doublereal>& inFrcbuf = pDM->GetBufIn(uInFrcLabel);
		const std::vector<doublereal>& outbuf = pDM->GetBufOut(uOutLabel);
#elif BUFFER_WHICH & BUFFER_RAW
		doublereal *inKinbuf;
		doublereal *inFrcbuf;
		const doublereal *outbuf;
#if BUFFER_WHICH & BUFFER_MBDYN_OWNS_MEM
		inKinbuf = pDM->GetBufInRaw(uInKinLabel);
		inFrcbuf = pDM->GetBufInRaw(uInFrcLabel);
		outbuf = pDM->GetBufOutRaw(uOutLabel);
#else
		inKinbuf = new doublereal[iInKinSize];
		inFrcbuf = new doublereal[iInFrcSize];
		outbuf = new doublereal[iOutSize];
		pDM->SetBufInRaw(uInKinLabel, iInKinSize, inKinbuf);
		pDM->SetBufInRaw(uInFrcLabel, iInFrcSize, inFrcbuf);
		pDM->SetBufOutRaw(uOutLabel, iOutSize, outbuf);
#endif
#endif

		const char *sInKinName = "prescribed_displacement.drv";
		std::ifstream streamInKin(sInKinName);
		if (!streamInKin) {
			std::cerr << "unable to open file \"" << sInKinName << "\"" << std::endl;
			exit(1);
		}

		const char *sInFrcName = "prescribed_force.drv";
		std::ifstream streamInFrc(sInFrcName);
		if (!streamInFrc) {
			std::cerr << "unable to open file \"" << sInFrcName << "\"" << std::endl;
			exit(1);
		}

		// write to I buffer(s)
		for (int i = 0; i < iInKinSize; i++) {
			streamInKin >> inKinbuf[i];
		}

		for (int i = 0; i < iInFrcSize; i++) {
			streamInFrc >> inFrcbuf[i];
		}
		// end of buffer initialization



		if (pS->Start()) {
			// read from O buffer
			if (!sOutFileName) {
				(*out) << "out[] =";
			}
			for (size_t i = 0; i < iOutSize; i++) {
				(*out) << " " << outbuf[i];
			}
			(*out) << std::endl;

			for (;;) {
				// write to I buffer(s)
				for (int i = 0; i < iInKinSize; i++) {
					streamInKin >> inKinbuf[i];
				}

				for (int i = 0; i < iInFrcSize; i++) {
					streamInFrc >> inFrcbuf[i];
				}

				// additional damping (beware: using velocity at past time step)
				inFrcbuf[0] -= c*outbuf[1];

				bool b = pS->Advance();

				// read from O buffer
				if (!sOutFileName) {
					(*out) << "out[] =";
				}
				for (size_t i = 0; i < iOutSize; i++) {
					(*out) << " " << outbuf[i];
				}
				(*out) << std::endl;

				if (!b) {
					break;
				}
			}
		}

		// cleanup I/O buffers
#if BUFFER_WHICH == BUFFER_RAW
		delete[] inKinbuf;
		delete[] inFrcbuf;
		delete[] outbuf;
#endif
	}

	delete pS;

	return 0;
}

