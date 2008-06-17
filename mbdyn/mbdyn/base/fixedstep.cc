/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <fstream>

#include "dataman.h"
#include "filedrv.h"
#include "fixedstep.h"
#include "solver.h"

/* FixedStepFileDrive - begin */

FixedStepFileDrive::FixedStepFileDrive(unsigned int uL,
		const DriveHandler* pDH,
		const char* const sFileName,
		integer ins, integer ind,
		doublereal t0, doublereal dt, bool pz, Drive::Bailout bo)
: FileDrive(uL, pDH, sFileName, ind),
dT0(t0), dDT(dt), iNumSteps(ins), bPadZeroes(pz), boWhen(bo), pd(0), pvd(0)
{
	ASSERT(iNumDrives > 0);
	ASSERT(sFileName != NULL);
	ASSERT(dDT > 0.);

	std::ifstream in(sFileName);
	if (!in) {
		silent_cerr("can't open file \""
			<< sFileName << "\"" << std::endl);
		throw ErrGeneric();
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

	if (ins == -1) {
		std::streampos pos = in.tellg();

		for (ins = 0; !in.eof(); ins++) {
			char buf[1024];

			do {
				in.getline(buf, sizeof(buf));
			} while (strlen(buf) == STRLENOF(buf) && buf[STRLENOF(buf)] != '\n');
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
				throw ErrGeneric();
			}
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
		if (bPadZeroes) {
			for (int i = 1; i <= iNumDrives; i++) {
				pdVal[i] = 0.;
			}

		} else if (boWhen & Drive::BO_LOWER) {
			throw Solver::EndOfSimulation(0,
				"A fixed step file drive lower bound is halting	the simulation");
		} else {
			for (int i = 1; i <= iNumDrives; i++) {
				pdVal[i] = pvd[i][0];
			}
		}

	} else if (tt > dDT*(iNumSteps - 1)) {
		if (bPadZeroes) {
			for (int i = 1; i <= iNumDrives; i++) {
				pdVal[i] = 0.;
			}

		} else if (boWhen & Drive::BO_UPPER) {
			throw Solver::EndOfSimulation(0,
				"A fixed step file drive upper bound is halting	the simulation");
		} else {
			for (int i = 1; i <= iNumDrives; i++) {
				pdVal[i] = pvd[i][iNumSteps - 1];
			}
		}

	} else {
		integer j1 = integer(floor(tt/dDT));
		if (j1 == iNumSteps) {
			j1--;
		}
		integer j2 = j1 + 1;
		doublereal dt1 = dT0 + j1*dDT;
		doublereal dt2 = dt1 + dDT;

		for (int i = 1; i <= iNumDrives; i++) {
   			pdVal[i] = (pvd[i][j2]*(t - dt1) - pvd[i][j1]*(t - dt2))/dDT;
		}
	}
}

/* FixedStepFileDrive - end */


/* legge i drivers tipo fixed step file */

Drive*
ReadFixedStepFileDrive(DataManager* pDM,
		MBDynParser& HP,
		unsigned int uLabel)
{
	integer isteps = -1;
	if (!HP.IsKeyWord("count")) {
		isteps = HP.GetInt();
		if (isteps <= 0) {
			silent_cerr("FixedStepFileDrive(" << uLabel << "): "
				"invalid steps number " << isteps
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric();
		}
	}

	integer idrives = HP.GetInt();
	if (idrives <= 0) {
		silent_cerr("FixedStepFileDrive(" << uLabel << "): "
			"invalid channels number " << idrives
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric();
	}

	doublereal t0 = HP.GetReal();
	doublereal dt = HP.GetReal();
	if (dt <= 0) {
		silent_cerr("FixedStepFileDrive(" << uLabel << "): "
			"invalid time step " << dt
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric();
	}

	bool pz(true);
	Drive::Bailout bo(Drive::BO_NONE);

	if (HP.IsKeyWord("pad" "zeros") || HP.IsKeyWord("pad" "zeroes")) {
		if (HP.IsKeyWord("no")) {
			pz = false;

		} else if (!HP.IsKeyWord("yes")) {
			silent_cerr("unknown value for \"pad zeros\" "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric();
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
			throw ErrGeneric();
		}
	}

	const char* filename = HP.GetFileName();

	Drive* pDr = NULL;
	SAFENEWWITHCONSTRUCTOR(pDr,
			FixedStepFileDrive,
			FixedStepFileDrive(uLabel, pDM->pGetDrvHdl(),
				filename, isteps, idrives,
				t0, dt, pz, bo));

	return pDr;
} /* End of ReadFixedStepFileDrive */

