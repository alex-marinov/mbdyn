/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2008
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>
#include "extforce.h"

#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

/* ExtForce - begin */

/* Costruttore */
ExtForce::ExtForce(unsigned int uL,
	std::string& fin,
	bool bRemoveIn,
        std::string& fout,
	bool bNoClobberOut,
	int iSleepTime,
	int iCoupling,
	int iPrecision,
	flag fOut)
: Elem(uL, fOut), 
Force(uL, fOut), 
fin(fin.c_str()),
fout(fout.c_str()),
bRemoveIn(bRemoveIn),
bNoClobberOut(bNoClobberOut),
bFirstRes(false),
iSleepTime(iSleepTime),
iCoupling(iCoupling),
iCouplingCounter(0),
iPrecision(iPrecision)
{
	NO_OP;
}

ExtForce::~ExtForce(void)
{
	NO_OP;
}


void
ExtForce::Update(const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	/* If running tight coupling, send kinematics every iteration */
	/* NOTE: tight coupling may need relaxation */
	if (iCoupling && !((++iCouplingCounter)%iCoupling)) {
		Send();
	}
}
	
/*
 * Elaborazione stato interno dopo la convergenza
 */
void
ExtForce::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	/* After prediction, mark next residual as first */
	bFirstRes = true;
}

/*
 * Elaborazione stato interno dopo la convergenza
 */
void
ExtForce::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP)
{
	/* If not running tight coupling, send kinematics only at convergence */
	if (iCoupling != 1) {
		Send();
		if (bRemoveIn) {
			Unlink();
		}
	}
}

/*
 * Unlink input file when no longer required
 * used to inform companion software that a new input file can be written
 */
void
ExtForce::Unlink(void)
{
	if (unlink(fin.c_str()) != 0) {
		int save_errno = errno;

		switch (save_errno) {
		case ENOENT:
			break;

		default:
			silent_cerr("ExtForce(" << GetLabel() << "): "
				<< "unable to delete input file \"" << fin.c_str() 
				<< "\": " << strerror(save_errno) << std::endl);
			throw ErrGeneric();
		}
	}
}

/*
 * Send output to companion software
 */
void
ExtForce::Send(void)
{
	if (bNoClobberOut) {
		bool	bKeepGoing(true);

		for (int cnt = 0; bKeepGoing; cnt++) {
			struct stat	s;

			if (stat(fout.c_str(), &s) != 0) {
				int save_errno = errno;

				switch (save_errno) {
				case ENOENT:
					bKeepGoing = false;
					break;

				default:
					silent_cerr("ExtForce(" << GetLabel() << "): "
						"unable to stat output file \"" << fout.c_str() << "\": "
						<< strerror(save_errno) << std::endl);
					throw ErrGeneric();
				}

			} else {
#ifdef USE_SLEEP
				silent_cout("ExtForce(" << GetLabel() << "): "
					"output file \"" << fout.c_str() << "\" still present, "
					"try #" << cnt << "; "
					"sleeping " << iSleepTime << " s" << std::endl);
				sleep(iSleepTime);
#endif // USE_SLEEP
			}
		}
	}

	std::string tmpout(fout + ".tmp");
	std::ofstream outf(tmpout.c_str());

	if (!outf) {
		silent_cerr("ExtForce(" << GetLabel() << "): "
			"unable to open file \"" << fout.c_str() << "\"" << std::endl);
		throw ErrGeneric();
	}

	if (iPrecision != 0) {
		outf.precision(iPrecision);
	}
	outf.setf(std::ios::scientific);

	Send(outf);

	/* send */
	outf.close();
	rename(tmpout.c_str(), fout.c_str());
}

void
ExtForce::Recv(void)
{
	if ((iCoupling && !(iCouplingCounter%iCoupling)) || bFirstRes) {
	        std::ifstream inf(fin.c_str());

#ifdef USE_SLEEP
		for (int cnt = 0; !inf; cnt++) {
			silent_cout("ExtForce(" << GetLabel() << "): "
				"input file \"" << fin.c_str() << "\" missing, "
				"try #" << cnt << "; "
				"sleeping " << iSleepTime << " s" << std::endl); 
               
			sleep(iSleepTime);
			inf.clear();
			inf.open(fin.c_str());
		}
#endif // USE_SLEEP

		Recv(inf);

		if (bRemoveIn) {
			Unlink();
		}
	}

	bFirstRes = false;
}

void
ReadExtForce(DataManager* pDM, 
	MBDynParser& HP, 
	unsigned int uLabel,
	std::string& fin, bool& bUnlinkIn,
	std::string& fout, bool& bNoClobberOut,
	int& iSleepTime,
	int& iCoupling,
	int& iPrecision)
{
	const char	*s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForce(" << uLabel << "): unable to get input file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}
	fin = s;

	bUnlinkIn = false;
	if (HP.IsKeyWord("unlink")) {
		bUnlinkIn = true;
	}

	s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForce(" << uLabel << "): unable to get output file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}
	fout = s;

	bNoClobberOut = false;
	if (HP.IsKeyWord("no" "clobber")) {
		bNoClobberOut = true;
	}

	iSleepTime = 1;
	if (HP.IsKeyWord("sleep" "time")) {
		iSleepTime = HP.GetInt();
		if (iSleepTime <= 0 ) {
			silent_cerr("ExtForce(" << uLabel << "): "
				"invalid sleep time " << iSleepTime <<std::endl);
			throw ErrGeneric();
		}
	}

	iCoupling = 0;
	if (HP.IsKeyWord("coupling")) {
		if (HP.IsKeyWord("loose")) {
			iCoupling = 0;
		} else if (HP.IsKeyWord("tight")) {
			iCoupling = 1;
		} else {
			iCoupling = HP.GetInt();
			if (iCoupling < 0) {
				silent_cerr("ExtForce(" << uLabel << "): "
					"invalid coupling value "
					"\"" << iCoupling << "\""
					" at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric();
			}
		}
	}

	iPrecision = 0;
	if (HP.IsKeyWord("precision")) {
		if (!HP.IsKeyWord("default")) {
			iPrecision = HP.GetInt();
		}

		if (iPrecision < 0) {
			silent_cerr("ExtForce(" << uLabel << "): "
				"invalid precision value "
				"\"" << iPrecision << "\""
				" at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric();
		}
	}
}

