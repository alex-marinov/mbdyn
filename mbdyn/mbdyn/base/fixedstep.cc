/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#include <ac/fstream>

#include <dataman.h>
#include <filedrv.h>
#include <fixedstep.h>

/* FixedStepFileDrive - begin */

FixedStepFileDrive::FixedStepFileDrive(unsigned int uL,
		const DriveHandler* pDH,
		const char* const sFileName,
		integer is, integer nd,
		doublereal t0, doublereal dt)
: FileDrive(uL, pDH, sFileName, nd),
dT0(t0), dDT(dt), iNumSteps(is), pd(NULL), pvd(NULL)
{
	ASSERT(iNumSteps > 0);
	ASSERT(iNumDrives > 0);
	ASSERT(sFileName != NULL);
	ASSERT(dDT > 0.);

	std::ifstream in(sFileName);
	if (!in) {
		silent_cerr("can't open file \""
			<< sFileName << "\"" << std::endl);
		throw ErrGeneric();
	}

	SAFENEWARR(pd, doublereal, iNumDrives*iNumSteps);
	SAFENEWARR(pvd, doublereal*, iNumDrives+1);

	/* Attenzione: il primo puntatore e' vuoto
	 * (ne e' stato allocato uno in piu'),
	 * cosi' i drives possono essere numerati da 1 a n */
	for (integer i = iNumDrives; i-- > 0; ) {
		pvd[i+1] = pd+i*iNumSteps;
	}

	/*
	 * Mangia gli eventuali commenti iniziali
	 */
	char c = '\0';
	while (in.get(c), c == '#') {
		char buf[1024];

		do {
			in.getline(buf, sizeof(buf));
		} while (strlen(buf) == sizeof(buf) - 1 && buf[sizeof(buf) - 1] != '\n');
	}

	if (c != '#') {
		in.putback(c);
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
	integer j1 = integer(floor((t-dT0)/dDT));
	integer j2 = j1+1;

	if (j2 < 0 || j1 > iNumSteps) {
		for (int i = 1; i <= iNumDrives; i++) {
			pdVal[i] = 0.;
		}

	} else {
		doublereal dt1 = dT0+j1*dDT;
		doublereal dt2 = dt1+dDT;

		for (int i = 1; i <= iNumDrives; i++) {
   			pdVal[i] = (pvd[i][j2]*(t-dt1)-pvd[i][j1]*(t-dt2))/dDT;
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
	integer isteps = HP.GetInt();
	integer idrives = HP.GetInt();
	doublereal t0 = HP.GetReal();
	doublereal dt = HP.GetReal();
	const char* filename = HP.GetFileName();

	Drive* pDr = NULL;
	SAFENEWWITHCONSTRUCTOR(pDr,
			FixedStepFileDrive,
			FixedStepFileDrive(uLabel, pDM->pGetDrvHdl(),
				filename, isteps, idrives,
				t0, dt));

	return pDr;
} /* End of ReadFixedStepFileDrive */

