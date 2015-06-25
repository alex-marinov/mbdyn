/* $Header$ */
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

/* fixed step file driver */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <fstream>

#include "dataman.h"
#include "filedrv.h"
#include "fixedstep.h"
#include "solver.h"

/* FixedStepFileDrive - begin */

static const std::vector<doublereal> v0;
static const doublereal dFromFile = -std::numeric_limits<doublereal>::max();

FixedStepFileDrive::FixedStepFileDrive(unsigned int uL,
		const DriveHandler* pDH,
		const char* const sFileName,
		integer ins, integer ind,
		doublereal t0, doublereal dt,
		bool bl, bool pz, Drive::Bailout bo)
: FileDrive(uL, pDH, sFileName, ind, v0),
dT0(t0), dDT(dt), iNumSteps(ins),
bLinear(bl), bPadZeroes(pz), boWhen(bo), pd(0), pvd(0)
{
	ASSERT(iNumDrives > 0);
	ASSERT(sFileName != NULL);
	ASSERT(dDT > 0.);

	std::ifstream in(sFileName);
	if (!in) {
		silent_cerr("FixedStepFileDrive(" << uL << "): "
			"can't open file \"" << sFileName << "\""
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/*
	 * Mangia gli eventuali commenti iniziali
	 */
	char cTmp = '\0';
	while (in.get(cTmp), cTmp == '#') {
		char tmpbuf[1024];

		do {
			in.getline(tmpbuf, sizeof(tmpbuf));
			int idx = 0;
			while (isspace(tmpbuf[idx])) {
				idx++;
			}

			if (dT0 == dFromFile && strncasecmp(&tmpbuf[idx], "initial time:", STRLENOF("initial time:")) == 0) {
				double d;
				if (sscanf(&tmpbuf[idx + STRLENOF("initial time:")], "%le", &d) != 1) {
					silent_cerr("FixedStepFileDrive(" << uL << "): "
						"can't parse \"initial time\" line \"" << &tmpbuf[idx] << "\" "
						"from file \"" << sFileName << "\"" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				dT0 = d;

			} else if (dDT == dFromFile && strncasecmp(&tmpbuf[idx], "time step:", STRLENOF("time step:")) == 0) {
				double d;
				if (sscanf(&tmpbuf[idx + STRLENOF("time step:")], "%le", &d) != 1) {
					silent_cerr("FixedStepFileDrive(" << uL << "): "
						"can't parse \"time step\" line \"" << &tmpbuf[idx] << "\" "
						"from file \"" << sFileName << "\"" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				if (d <= 0.) {
					silent_cerr("FixedStepFileDrive(" << uL << "): "
						"invalid time step value " << d << " "
						"in file \"" << sFileName << "\"" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				dDT = d;
			}
		} while (strlen(tmpbuf) == STRLENOF(tmpbuf) && tmpbuf[STRLENOF(tmpbuf)] != '\n');
	}

	if (dT0 == dFromFile) {
		silent_cerr("FixedStepFileDrive(" << uL << "): "
			"expecting \"initial time\" line in file \"" << sFileName << "\""
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (dDT == dFromFile) {
		silent_cerr("FixedStepFileDrive(" << uL << "): "
			"expecting \"time step\" line in file \"" << sFileName << "\""
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (cTmp != '#') {
		in.putback(cTmp);
	}

	if (ins == -1) {
		std::streampos pos = in.tellg();

		for (ins = 0; !in.eof(); ins++) {
			char tmpbuf[1024];

			do {
				in.getline(tmpbuf, sizeof(tmpbuf));
			} while (strlen(tmpbuf) == STRLENOF(tmpbuf) && tmpbuf[STRLENOF(tmpbuf)] != '\n');
		}
		iNumSteps = --ins;

		in.clear();
		in.seekg(pos);

		silent_cout("FixedStepFileDrive(" << uL << "): "
			"counted " << ins << " steps" << std::endl);
	}

	SAFENEWARR(pd, doublereal, iNumDrives*iNumSteps);
	SAFENEWARR(pvd, doublereal*, iNumDrives + 1);

	/* Attenzione: il primo puntatore e' vuoto
	 * (ne e' stato allocato uno in piu'),
	 * cosi' i drives possono essere numerati da 1 a n */
	for (integer i = iNumDrives; i-- > 0; ) {
		pvd[i + 1] = pd + i*iNumSteps;
	}

	for (integer j = 0; j < iNumSteps; j++) {
		for (integer i = 1; i <= iNumDrives; i++) {
			in >> pvd[i][j];
			if (in.eof()) {
				silent_cerr("unexpected end of file '"
					<< sFileName << '\'' << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			int c;
			for (c = in.get(); isspace(c); c = in.get()) {
				if (c == '\n') {
					if (i != iNumDrives) {
						silent_cerr("unexpected end of line #" << j + 1 << " after channel #" << i << ", column #" << i << " of file '"
							<< sFileName << '\'' << std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
					break;
				}
			}

			if (i == iNumDrives && c != '\n') {
				silent_cerr("missing end-of-line at line #" << j + 1 << " of file '"
					<< sFileName << '\'' << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			in.putback(c);
		}
	}

	/* All data is available, so initialize the buffer accordingly */
	ServePending(pDH->dGetTime());
}

FixedStepFileDrive::~FixedStepFileDrive(void)
{
	SAFEDELETEARR(pd);
	SAFEDELETEARR(pvd);
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
FixedStepFileDrive::Restart(std::ostream& out) const
{
	return out << "0. /* FixedStepFileDrive: not implemented yet! */"
		<< std::endl;
}

void
FixedStepFileDrive::ServePending(const doublereal& t)
{
	doublereal tt = t - dT0;

	if (tt < 0) {
		if (boWhen & Drive::BO_LOWER) {
			throw Solver::EndOfSimulation(EXIT_SUCCESS,
				MBDYN_EXCEPT_ARGS,
				"A fixed step file drive lower bound is halting the simulation");
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

	} else if (tt > dDT*(iNumSteps - 1)) {
		if (boWhen & Drive::BO_UPPER) {
			throw Solver::EndOfSimulation(EXIT_SUCCESS,
				MBDYN_EXCEPT_ARGS,
				"A fixed step file drive upper bound is halting the simulation");
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
		integer j1 = integer(floor(tt/dDT));
		if (bLinear) {
			if (j1 == iNumSteps) {
				j1--;
			}
			integer j2 = j1 + 1;
			doublereal dt1 = dT0 + j1*dDT;
			doublereal dt2 = dt1 + dDT;
			doublereal dw1 = (dt2 - t)/dDT;
			doublereal dw2 = (t - dt1)/dDT;

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

/* FixedStepFileDrive - end */


/* legge i drivers tipo fixed step file */

Drive *
FixedStepDR::Read(unsigned uLabel, const DataManager *pDM, MBDynParser& HP)
{
	integer isteps = -1;
	if (!HP.IsKeyWord("count")) {
		isteps = HP.GetInt();
		if (isteps <= 0) {
			silent_cerr("FixedStepFileDrive(" << uLabel << "): "
				"invalid steps number " << isteps
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	integer idrives = HP.GetInt();
	if (idrives <= 0) {
		silent_cerr("FixedStepFileDrive(" << uLabel << "): "
			"invalid channels number " << idrives
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!HP.IsKeyWord("initial" "time")) {
		silent_cerr("FixedStepFileDrive(" << uLabel << "): "
			"warning, keyword \"initial time\" expected at line "
		<< HP.GetLineData() << std::endl);
	}
	doublereal t0 = dFromFile;
	if (!HP.IsKeyWord("from" "file")) {
		t0 = HP.GetReal();
	}

	if (!HP.IsKeyWord("time" "step")) {
		silent_cerr("FixedStepFileDrive(" << uLabel << "): "
			"warning, keyword \"time step\" expected at line "
		<< HP.GetLineData() << std::endl);
	}
	doublereal dt = dFromFile;
	if (!HP.IsKeyWord("from" "file")) {
		dt = HP.GetReal();
		if (dt <= 0.) {
			silent_cerr("FixedStepFileDrive(" << uLabel << "): "
				"invalid time step " << dt
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	bool bl(true);
	if (HP.IsKeyWord("interpolation")) {
		if (HP.IsKeyWord("const")) {
			bl = false;

		} else if (!HP.IsKeyWord("linear")) {
			silent_cerr("FixedStepFileDrive(" << uLabel << "): "
				"unknown value for \"interpolation\" "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	bool pz(true);
	Drive::Bailout bo(Drive::BO_NONE);

	if (HP.IsKeyWord("pad" "zeros") || HP.IsKeyWord("pad" "zeroes")) {
		if (!HP.GetYesNo(pz)) {
			silent_cerr("FixedStepFileDrive(" << uLabel << "): "
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
			silent_cerr("FixedStepFileDrive(" << uLabel << "): "
				"invalid bailout parameter "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	const char* filename = HP.GetFileName();

	Drive* pDr = NULL;
	SAFENEWWITHCONSTRUCTOR(pDr,
			FixedStepFileDrive,
			FixedStepFileDrive(uLabel, pDM->pGetDrvHdl(),
				filename, isteps, idrives,
				t0, dt, bl, pz, bo));

	return pDr;
}

