/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/* fixed step file driver */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <fstream>

#include "dataman.h"
#include "filedrv.h"
#include "varstep.h"
#include "solver.h"
#include "bisec.h"

/* VariableStepFileDrive - begin */

static const std::vector<doublereal> v0;

VariableStepFileDrive::VariableStepFileDrive(unsigned int uL,
		const DriveHandler* pDH,
		const char* const sFileName,
		integer ind, bool bl, bool pz, Drive::Bailout bo)
: FileDrive(uL, pDH, sFileName, ind, v0),
iNumSteps(-1), iCurrStep(-1),
bLinear(bl), bPadZeroes(pz), boWhen(bo), pd(0), pvd(0)
{
	ASSERT(iNumDrives > 0);
	ASSERT(sFileName != NULL);

	std::ifstream in(sFileName);
	if (!in) {
		silent_cerr("can't open file \""
			<< sFileName << "\"" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/*
	 * Mangia gli eventuali commenti iniziali
	 */
	char c = '\0';
	while (in.get(c), c == '#') {
		char buf[1024];

		do {
			in.getline(buf, sizeof(buf));
		} while (strlen(buf) == STRLENOF(buf) && buf[STRLENOF(buf)] != '\n');
	}

	if (c != '#') {
		in.putback(c);
	}

	if (iNumSteps == -1) {
		std::streampos pos = in.tellg();

		integer ins;
		for (ins = 0; !in.eof(); ins++) {
			char buf[1024];

			do {
				in.getline(buf, sizeof(buf));
			} while (strlen(buf) == STRLENOF(buf) && buf[STRLENOF(buf)] != '\n');
		}
		iNumSteps = --ins;

		in.clear();
		in.seekg(pos);

		silent_cout("VariableStepFileDrive(" << uL << "): "
			"counted " << ins << " steps" << std::endl);
	}

	SAFENEWARR(pd, doublereal, (1 + iNumDrives)*iNumSteps);
	SAFENEWARR(pvd, doublereal*, 1 + iNumDrives);

	for (integer i = iNumDrives + 1; i-- > 0; ) {
		pvd[i] = pd + i*iNumSteps;
	}

	for (integer j = 0; j < iNumSteps; j++) {
		// 0 -> iNumDrives to account for time
		for (integer i = 0; i <= iNumDrives; i++) {
			in >> pvd[i][j];
			if (in.eof()) {
				silent_cerr("unexpected end of file '"
					<< sFileName << '\'' << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			int c;
			for (c = in.get(); isspace(c); c = in.get()) {
				if (c == '\n') {
					if (i == 0) {
						silent_cerr("unexpected end of line #" << j + 1 << " after time=" << pvd[0][j] << ", column #" << i + 1 << " of file '"
							<< sFileName << '\'' << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					if (i != iNumDrives) {
						silent_cerr("unexpected end of line #" << j + 1 << ", time=" << pvd[0][j] << " after channel #" << i << ", column #" << i + 1 << " of file '"
							<< sFileName << '\'' << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					break;
				}
			}

			if (i == iNumDrives && c != '\n') {
				silent_cerr("missing end-of-line at line #" << j + 1 << ", time=" << pvd[0][j] << " of file '"
					<< sFileName << '\'' << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			in.putback(c);
		}

		if (j > 0) {
			if (pvd[0][j] <= pvd[0][j - 1]) {
				silent_cerr("time[" << j << "]=" << pvd[0][j]
					<< " <= time[" << j - 1 << "]=" << pvd[0][j - 1]
					<< " in file '" << sFileName << "'" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	// All data is available, so initialize the buffer accordingly
	// use bisection to initialize iCurrStep
	doublereal dTime = pDH->dGetTime();
	iCurrStep = bisec(pvd[0], dTime, 0, iNumSteps - 1);
	if (iCurrStep < 0) {
		iCurrStep++;
	}

	ServePending(dTime);
}

VariableStepFileDrive::~VariableStepFileDrive(void)
{
	SAFEDELETEARR(pd);
	SAFEDELETEARR(pvd);
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
VariableStepFileDrive::Restart(std::ostream& out) const
{
	return out << "0. /* VariableStepFileDrive: not implemented yet! */"
		<< std::endl;
}

void
VariableStepFileDrive::ServePending(const doublereal& t)
{
	if (t <= pvd[0][0]) {
		if (t < pvd[0][0] && (boWhen & Drive::BO_LOWER)) {
			throw Solver::EndOfSimulation(EXIT_SUCCESS,
				MBDYN_EXCEPT_ARGS,
				"A variable step file drive lower bound is halting the simulation");
		}

		if (bPadZeroes) {
			for (int i = 1; i <= iNumDrives; i++) {
				pdVal[i] = 0.;
			}

		} else {
			for (int i = 1; i <= iNumDrives; i++) {
				pdVal[i] = pvd[i][0];
			}
		}

	} else if (t >= pvd[0][iNumSteps - 1]) {
		if (t > pvd[0][iNumSteps - 1] && (boWhen & Drive::BO_UPPER)) {
			throw Solver::EndOfSimulation(EXIT_SUCCESS,
				MBDYN_EXCEPT_ARGS,
				"A variable step file drive upper bound is halting the simulation");
		}

		if (bPadZeroes) {
			for (int i = 1; i <= iNumDrives; i++) {
				pdVal[i] = 0.;
			}

		} else {
			for (int i = 1; i <= iNumDrives; i++) {
				pdVal[i] = pvd[i][iNumSteps - 1];
			}
		}

	} else {
		// look for step exactly before
		// note: use linear search, under the assumption
		// that we always start from a relatively close guess
		while (pvd[0][iCurrStep] > t) {
			iCurrStep--;
		}

		// no need to check for under/uverflow
		ASSERT(iCurrStep >= 0);

		while (pvd[0][iCurrStep + 1] <= t) {
			iCurrStep++;
		}
		ASSERT(iCurrStep < iNumSteps);

		integer j1 = iCurrStep;
		if (bLinear) {
			integer j2 = j1 + 1;
			doublereal dt1 = pvd[0][j1];
			doublereal dt2 = pvd[0][j2];
			doublereal dw1 = (dt2 - t)/(dt2 - dt1);
			doublereal dw2 = (t - dt1)/(dt2 - dt1);

			for (int i = 1; i <= iNumDrives; i++) {
   				pdVal[i] = pvd[i][j2]*dw2 + pvd[i][j1]*dw1;
			}

		} else {
			for (int i = 1; i <= iNumDrives; i++) {
   				pdVal[i] = pvd[i][j1];
			}
		}
	}
}

/* VariableStepFileDrive - end */


/* legge i drivers tipo fixed step file */

Drive *
VariableStepDR::Read(unsigned uLabel, const DataManager *pDM, MBDynParser& HP)
{
	integer idrives = HP.GetInt();
	if (idrives <= 0) {
		silent_cerr("VariableStepFileDrive(" << uLabel << "): "
			"invalid channels number " << idrives
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	bool bl(true);
	if (HP.IsKeyWord("interpolation")) {
		if (HP.IsKeyWord("const")) {
			bl = false;

		} else if (!HP.IsKeyWord("linear")) {
			silent_cerr("VariableStepFileDrive(" << uLabel << "): "
				"unknown value for \"interpolation\" "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	bool pz(true);
	Drive::Bailout bo(Drive::BO_NONE);

	if (HP.IsKeyWord("pad" "zeros") || HP.IsKeyWord("pad" "zeroes")) {
		if (!HP.GetYesNo(pz)) {
			silent_cerr("VariableStepFileDrive(" << uLabel << "): "
				"unknown value for \"pad zeros\" "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else if (HP.IsKeyWord("bailout")) {
		if (HP.IsKeyWord("none")) {
			bo = Drive::BO_NONE;

		} else if (HP.IsKeyWord("upper")) {
			bo = Drive::BO_UPPER;

		} else if (HP.IsKeyWord("lower")) {
			bo = Drive::BO_LOWER;

		} else if (HP.IsKeyWord("any")) {
			bo = Drive::BO_ANY;

		} else {
			silent_cerr("VariableStepFileDrive(" << uLabel << "): "
				"invalid bailout parameter "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	const char* filename = HP.GetFileName();

	Drive* pDr = NULL;
	SAFENEWWITHCONSTRUCTOR(pDr,
			VariableStepFileDrive,
			VariableStepFileDrive(uLabel, pDM->pGetDrvHdl(),
				filename, idrives, bl, pz, bo));

	return pDr;
}

