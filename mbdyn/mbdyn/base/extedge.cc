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

#include "dataman.h"
#include "extedge.h"
#include "except.h"
#include "solver.h"

#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

/* ExtFileHandlerEDGE - begin */

ExtFileHandlerEDGE::ExtFileHandlerEDGE(std::string& fflagname,
	std::string& fdataname, int iSleepTime, int iPrecision)
: ExtFileHandlerBase(iSleepTime, iPrecision),
fflagname(fflagname), fdataname(fdataname),
bReadForces(true)
{
	NO_OP;
}

ExtFileHandlerEDGE::~ExtFileHandlerEDGE(void)
{
	NO_OP;
}

ExtFileHandlerEDGE::EDGEcmd
ExtFileHandlerEDGE::CheckFlag(int& cnt)
{
	int cmd;

	cnt++;

	infile.open(fflagname.c_str());

#ifdef USE_SLEEP
	if (iSleepTime > 0) {
		// WARNING: loops forever
		// add optional, configurable limit?
		for (; !infile; cnt++) {
			silent_cout("flag file \"" << fflagname.c_str() << "\" "
				"missing, try #" << cnt << "; "
				"sleeping " << iSleepTime << " s" << std::endl); 
			if (mbdyn_stop_at_end_of_iteration()) {
				cmd = EDGE_QUIT;
				goto done;
			}

			sleep(iSleepTime);

			infile.clear();
			infile.open(fflagname.c_str());
		}
	} else
#endif // USE_SLEEP
	{
		// error
		silent_cerr("flag file \"" << fflagname.c_str() << "\" "
			"missing" << std::endl); 
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	infile >> cmd;

done:;
	silent_cerr("flag file \"" << fflagname.c_str() << "\": cmd=" << cmd << std::endl);

	infile.close();
	infile.clear();

	return EDGEcmd(cmd);
}

void
ExtFileHandlerEDGE::AfterPredict(void)
{
	bReadForces = true;
}

std::ostream&
ExtFileHandlerEDGE::Send_pre(bool bAfterConvergence)
{
	int cnt = 0;
	if (!bReadForces) {
		goto bad;
	}

	// open data file for writing
	switch (CheckFlag(cnt)) {
	case EDGE_READ_READY:
	case EDGE_GOTO_NEXT_STEP:
		outfile.open(fdataname.c_str());
		if (!outfile) {
			silent_cerr("unable to open data file "
				"\"" << fdataname.c_str() << "\" "
				"for output" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		outfile.setf(std::ios::scientific);
		outfile.precision(iPrecision);
		break;

	default:
bad:;
		outfile.setstate(std::ios_base::badbit);
		break;
	}

	return outfile;
}

void
ExtFileHandlerEDGE::Send_post(bool bAfterConvergence)
{
	if (!bReadForces) {
		return;
	}

	// close data file after writing
	if (outfile.is_open()) {
		outfile.close();
	}

	if (bAfterConvergence) {
#if 0
		// This stops EDGE's subiterations
		// and advances to next step
		outfile.open("update.ainp");
		outfile <<
			"UPDATE,N,0,0,1\n"
			"IBREAK,I,1,1,0\n"
			"5\n";
		outfile.close();
#endif
	}

	// write to a temporary, unique file and then rename to fflagname
#ifdef HAVE_MKSTEMP
	char ftmpname[] = "mbedgeXXXXXX";
	int filedes = mkstemp(ftmpname);
	FILE *fd = fdopen(filedes, "w");
#else // ! HAVE_MKSTEMP
	std::string ftmpn(fflagname + ".tmp");
	const char *ftmpname = ftmpn.c_str();
	FILE *fd = fopen(ftmpname, "w");
#endif // ! HAVE_MKSTEMP

	fprintf(fd, "%d", int(EDGE_MBDYN_WRITE_DONE));
	fclose(fd);
retry:;
	if (rename(ftmpname, fflagname.c_str()) == -1) {
		switch(errno) {
		case EBUSY:
#ifdef USE_SLEEP
			// TODO: configurable?
			sleep(1);
#endif // USE_SLEEP
			goto retry;
		default:
			silent_cerr("unable to rename flag file "
				"\"" << fdataname.c_str() << "\" "
				"(errno=" << errno << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	outfile.open(fflagname.c_str());
	outfile << int(EDGE_MBDYN_WRITE_DONE);
	outfile.close();
	outfile.clear();
}

std::istream&
ExtFileHandlerEDGE::Recv_pre(void)
{
	for (int cnt = 0; ;) {

		EDGEcmd cmd = CheckFlag(cnt);
		switch (cmd) {
		case EDGE_READ_READY:
			goto done;

		case EDGE_GOTO_NEXT_STEP:
			// FIXME: EDGE is done; do not read forces,
			// keep using old
			bReadForces = false;
			goto done;

		case EDGE_QUIT:
			silent_cout("EDGE requested end of simulation"
				<< std::endl);
			throw NoErr(MBDYN_EXCEPT_ARGS);

		default:
			break;
		}

		if (mbdyn_stop_at_end_of_iteration()) {
			bReadForces = false;
			goto done;
		}

		// WARNING: loops forever
		// add optional, configurable limit?
#ifdef USE_SLEEP
		if (iSleepTime > 0) {
			silent_cout("flag file \"" << fflagname.c_str() << "\": "
				"cmd=" << cmd << " try #" << cnt << "; "
				"sleeping " << iSleepTime << " s" << std::endl); 

			sleep(iSleepTime);
		} else
#endif // USE_SLEEP
		{
			silent_cout("flag file \"" << fflagname.c_str() << "\": "
				"cmd=" << cmd << " try #" << cnt << std::endl); 
		}
	}

done:;
	if (bReadForces) {
		infile.open(fdataname.c_str());
		if (!infile) {
			silent_cerr("unable to open data file "
				"\"" << fdataname.c_str() << "\" "
				"for input" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		infile.setstate(std::ios_base::badbit);
	}

	return infile;
}

bool
ExtFileHandlerEDGE::Recv_post(void)
{
	if (infile.is_open()) {
		infile.close();
	}

	return !bReadForces;
}

/* ExtFileHandlerEDGE - end */

ExtFileHandlerBase *
ReadExtFileHandlerEDGE(DataManager* pDM,
	MBDynParser& HP, 
	unsigned int uLabel)
{
	ExtFileHandlerBase *pEFH = 0;

	const char *s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForceEDGE(" << uLabel << "): "
			"unable to get flag file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	std::string fflagname = s;

	s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForceEDGE(" << uLabel << "): "
			"unable to get data file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	std::string fdataname = s;

	int iSleepTime = 1;
	int iPrecision = 0;
	ReadExtFileParams(pDM, HP, uLabel, iSleepTime, iPrecision);

	SAFENEWWITHCONSTRUCTOR(pEFH, ExtFileHandlerEDGE,
		ExtFileHandlerEDGE(fflagname, fdataname, iSleepTime, iPrecision));

	return pEFH;
}

