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
#include "except.h"

#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

/* ExtFileHandlerBase - begin */

ExtFileHandlerBase::~ExtFileHandlerBase(void)
{
	NO_OP;
}

/* ExtFileHandlerBase - end */

/* ExtFileHandler - begin */

ExtFileHandler::ExtFileHandler(std::string& fin,
	bool bRemoveIn,
        std::string& fout,
	bool bNoClobberOut,
	int iSleepTime,
	int iPrecision)
: fin(fin), fout(fout), tmpout(fout + ".tmp"),
bRemoveIn(bRemoveIn), bNoClobberOut(bNoClobberOut),
iSleepTime(iSleepTime), iPrecision(iPrecision)
{
	NO_OP;
}

ExtFileHandler::~ExtFileHandler(void)
{
	NO_OP;
}

std::ostream&
ExtFileHandler::Send_pre(bool bAfterConvergence)
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
					silent_cerr("unable to stat "
						"output file "
						"\"" << fout.c_str() << "\": "
						<< strerror(save_errno)
						<< std::endl);
					throw ErrGeneric();
				}

			} else {
#ifdef USE_SLEEP
				silent_cout("output file "
					"\"" << fout.c_str() << "\" "
					"still present, "
					"try #" << cnt << "; "
					"sleeping " << iSleepTime << " s"
					<< std::endl);
				sleep(iSleepTime);
#endif // USE_SLEEP
			}
		}
	}

	outf.open(tmpout.c_str());

	if (!outf) {
		silent_cerr("unable to open file \"" << fout.c_str() << "\""
			<< std::endl);
		throw ErrGeneric();
	}

	if (iPrecision != 0) {
		outf.precision(iPrecision);
	}
	outf.setf(std::ios::scientific);
	return outf;
}

void
ExtFileHandler::Send_post(bool bAfterConvergence)
{
	outf.close();
	rename(tmpout.c_str(), fout.c_str());
}

std::istream&
ExtFileHandler::Recv_pre(void)
{
	inf.open(fin.c_str());

#ifdef USE_SLEEP
	for (int cnt = 0; !inf; cnt++) {
		silent_cout("input file \"" << fin.c_str() << "\" missing, "
			"try #" << cnt << "; "
			"sleeping " << iSleepTime << " s" << std::endl); 
               
		sleep(iSleepTime);
		inf.clear();
		inf.open(fin.c_str());
	}
#endif // USE_SLEEP

	return inf;
}

void
ExtFileHandler::Recv_post(void)
{
	if (bRemoveIn) {
		if (unlink(fin.c_str()) != 0) {
			int save_errno = errno;

			switch (save_errno) {
			case ENOENT:
				break;

			default:
				silent_cerr("unable to delete input file "
					"\"" << fin.c_str() << "\": "
					<< strerror(save_errno) << std::endl);
				throw ErrGeneric();
			}
		}
	}
}

/* ExtFileHandler - end */

/* ExtFileHandlerEDGE - begin */

ExtFileHandlerEDGE::ExtFileHandlerEDGE(std::string& fflagname,
	std::string& fdataname, int iSleepTime)
: fflagname(fflagname), fdataname(fdataname), iSleepTime(iSleepTime)
{
	NO_OP;
}

ExtFileHandlerEDGE::~ExtFileHandlerEDGE(void)
{
	NO_OP;
}

std::ostream&
ExtFileHandlerEDGE::Send_pre(bool bAfterConvergence)
{
	outf.open(fdataname.c_str());
	if (!outf) {
		silent_cerr("unable to open data file "
			"\"" << fdataname.c_str() << "\" "
			"for output" << std::endl);
		throw ErrGeneric();
	}

	outf.setf(std::ios::scientific);
	return outf;
}

void
ExtFileHandlerEDGE::Send_post(bool bAfterConvergence)
{
	outf.close();

	if (bAfterConvergence) {
		outf.open("update.ainp");
		outf << "UPDATE,N,0,0,1\n"
			"IBREAK,I,1,1,0\n"
			"5\n";

	} else {
		outf.open(fflagname.c_str());
		outf << 3 << std::endl;
	}
	outf.close();
}

std::istream&
ExtFileHandlerEDGE::Recv_pre(void)
{
	for (int cnt = 0; ; cnt++) {
		inf.open(fflagname.c_str());

#ifdef USE_SLEEP
		for (; !inf; cnt++) {
			silent_cout("flag file \"" << fflagname.c_str() << "\" "
				"missing, try #" << cnt << "; "
				"sleeping " << iSleepTime << " s" << std::endl); 
               
			sleep(iSleepTime);
			inf.clear();
			inf.open(fdataname.c_str());
		}
#endif // USE_SLEEP

		int cmd;
		inf >> cmd;

		inf.close();

		if (cmd == 2) {
			break;

		} else if (cmd == 4) {
			silent_cout("EDGE requested end of simulation"
				<< std::endl);
			throw NoErr();
		}

		silent_cout("flag file \"" << fflagname.c_str() << "\": "
			"cmd=" << cmd << " try #" << cnt << "; "
			"sleeping " << iSleepTime << " s" << std::endl); 
#ifdef USE_SLEEP
		sleep(iSleepTime);
#endif // USE_SLEEP
	}

	inf.open(fdataname.c_str());
	if (!inf) {
		silent_cerr("unable to open data file "
			"\"" << fdataname.c_str() << "\" "
			"for input" << std::endl);
		throw ErrGeneric();
	}

	return inf;
}

void
ExtFileHandlerEDGE::Recv_post(void)
{
	inf.close();
}

/* ExtFileHandlerEDGE - end */

/* ExtForce - begin */

/* Costruttore */
ExtForce::ExtForce(unsigned int uL,
	ExtFileHandlerBase *pEFH,
	int iCoupling,
	flag fOut)
: Elem(uL, fOut), 
Force(uL, fOut),
pEFH(pEFH),
bFirstRes(false),
iCoupling(iCoupling),
iCouplingCounter(0)
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
		Send(true);
#if 0
		if (bRemoveIn) {
			Unlink();
		}
#endif
	}
}

/*
 * Send output to companion software
 */
void
ExtForce::Send(bool bAfterConvergence)
{
	std::ostream& outf = pEFH->Send_pre(bAfterConvergence);
	Send(outf, bAfterConvergence);
	pEFH->Send_post(bAfterConvergence);
}

void
ExtForce::Recv(void)
{
	if ((iCoupling && !(iCouplingCounter%iCoupling)) || bFirstRes) {
		std::istream& inf = pEFH->Recv_pre();
		Recv(inf);
		pEFH->Recv_post();
	}

	bFirstRes = false;
}

void
ExtForce::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0; 
	*piNumCols = 0; 
}
   
/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ExtForce::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ExtForce::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	WorkVec.ResizeReset(0);
	return WorkVec;
}

static ExtFileHandlerBase *
ReadExtFileHandler(DataManager* pDM,
	MBDynParser& HP, 
	unsigned int uLabel)
{
	ExtFileHandlerBase *pEFH = 0;

	if (HP.IsKeyWord("EDGE")) {
		
		return pEFH;
	}

	// default
	const char	*s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForce(" << uLabel << "): "
			"unable to get input file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}
	std::string fin = s;

	bool bUnlinkIn = false;
	if (HP.IsKeyWord("unlink")) {
		bUnlinkIn = true;
	}

	s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForce(" << uLabel << "): "
			"unable to get output file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}
	std::string fout = s;

	bool bNoClobberOut = false;
	if (HP.IsKeyWord("no" "clobber")) {
		bNoClobberOut = true;
	}

	int iSleepTime = 1;
	if (HP.IsKeyWord("sleep" "time")) {
		iSleepTime = HP.GetInt();
		if (iSleepTime <= 0 ) {
			silent_cerr("ExtForce(" << uLabel << "): "
				"invalid sleep time " << iSleepTime
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric();
		}
	}

	int iPrecision = 0;
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

	SAFENEWWITHCONSTRUCTOR(pEFH, ExtFileHandler,
		ExtFileHandler(fin, bUnlinkIn, fout, bNoClobberOut,
			iSleepTime, iPrecision));

	return pEFH;
}

void
ReadExtForce(DataManager* pDM, 
	MBDynParser& HP, 
	unsigned int uLabel,
	ExtFileHandlerBase*& pEFH,
	int& iCoupling)
{
	pEFH = ReadExtFileHandler(pDM, HP, uLabel);

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
}

