/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2017
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "extforce.h"
#include "extedge.h"
#include "extsocket.h"
#include "extsharedmem.h"
#include "except.h"
#include "solver.h"

#include <fstream>
#include <cstdlib>
#include <cerrno>
#include <sys/stat.h>

/* ExtFileHandlerBase - begin */

ExtFileHandlerBase::ExtFileHandlerBase(mbsleep_t SleepTime,
	std::streamsize Precision)
: Precision(Precision), SleepTime(SleepTime), bOK(true)
{
	NO_OP;
}

ExtFileHandlerBase::~ExtFileHandlerBase(void)
{
	NO_OP;
}

ExtFileHandlerBase::Negotiate
ExtFileHandlerBase::NegotiateRequest(void) const
{
	return NEGOTIATE_NO;
}

/* NOTE: getting here, in general, should be considered Bad (TM)
 * however, right now, it is used to distinguish whether communication
 * will occur on a iostream or a file descriptor */
std::ostream *
ExtFileHandlerBase::GetOutStream(void)
{
	return 0;
}

std::istream *
ExtFileHandlerBase::GetInStream(void)
{
	return 0;
}

#ifdef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP
mbdyn::shared_memory_buffer*
ExtFileHandlerBase::GetSharedMemBuffer(void)
{
	return 0;
}
#endif // HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

/* NOTE: getting here Is Bad (TM) */
int
ExtFileHandlerBase::GetOutFileDes(void)
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	// return -1;
}

int
ExtFileHandlerBase::GetSendFlags(void) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	// return 0;
}

int
ExtFileHandlerBase::GetInFileDes(void)
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	// return -1;
}

int
ExtFileHandlerBase::GetRecvFlags(void) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	// return 0;
}

/* ExtFileHandlerBase - end */

/* ExtFileHandler - begin */

ExtFileHandler::ExtFileHandler(std::string& fin,
	bool bRemoveIn,
        std::string& fout,
	bool bNoClobberOut,
	mbsleep_t SleepTime,
	std::streamsize Precision)
: ExtFileHandlerBase(SleepTime, Precision),
fin(fin), fout(fout), tmpout(fout + ".tmp"),
bRemoveIn(bRemoveIn), bNoClobberOut(bNoClobberOut)
{
	NO_OP;
}

ExtFileHandler::~ExtFileHandler(void)
{
	NO_OP;
}

bool
ExtFileHandler::Prepare_pre(void)
{
	return true;
}

void
ExtFileHandler::Prepare_post(bool ok)
{
	NO_OP;
}

void
ExtFileHandler::AfterPredict(void)
{
	NO_OP;
}

bool
ExtFileHandler::Send_pre(SendWhen when)
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
					silent_cerr("unable to stat output file "
						"\"" << fout.c_str() << "\": "
						<< strerror(save_errno)
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

			} else {
				if (mbdyn_stop_at_end_of_iteration()) {
					outf.setstate(std::ios_base::badbit);
					return false;
				}

				if (SleepTime > 0) {
					silent_cout("output file "
						"\"" << fout.c_str() << "\" "
						"still present, "
						"try #" << cnt << "; "
						"sleeping " << SleepTime << " s"
						<< std::endl);
					mbsleep(&SleepTime);

				} else {
					silent_cout("output file "
						"\"" << fout.c_str() << "\" "
						"still present, "
						"try #" << cnt
						<< std::endl);
				}
			}
		}
	}

	outf.open(tmpout.c_str());

	if (!outf) {
		silent_cerr("unable to open file \"" << fout.c_str() << "\""
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (Precision != 0) {
		outf.precision(Precision);
	}
	outf.setf(std::ios::scientific);
	return outf.good();
}

void
ExtFileHandler::Send_post(SendWhen when)
{
	outf.close();
	if (rename(tmpout.c_str(), fout.c_str()) != 0) {
		int save_errno = WSAGetLastError();
		silent_cerr("ExtFileHandler: unable to rename output file "
			"\"" << tmpout.c_str() << "\" "
			"into \"" << fout.c_str() << "\" "
			"(" << save_errno << ": " << strerror(save_errno) << ")"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

bool
ExtFileHandler::Recv_pre(void)
{
	inf.open(fin.c_str());

	for (int cnt = 0; !inf; cnt++) {
		silent_cout("input file \"" << fin.c_str() << "\" missing, "
			"try #" << cnt << "; "
			"sleeping " << SleepTime << " s" << std::endl); 
               
		if (mbdyn_stop_at_end_of_iteration()) {
			inf.setstate(std::ios_base::badbit);
			return (bOK = false);
		}

		mbsleep(&SleepTime);
		inf.clear();
		inf.open(fin.c_str());
	}

	return inf.good();
}

bool
ExtFileHandler::Recv_post(void)
{
	inf.close();

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
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	// NOTE: allow MBDyn to decide when converged
	return true;
}

std::ostream *
ExtFileHandler::GetOutStream(void)
{
	return &outf;
}

std::istream *
ExtFileHandler::GetInStream(void)
{
	return &inf;
}

/* ExtFileHandler - end */

/* ExtRemoteHandler - begin */

ESCmd
ExtRemoteHandler::u2cmd(unsigned u) const
{
	switch (u) {
	case ES_REGULAR_DATA:
	case ES_GOTO_NEXT_STEP:
	case ES_ABORT:
	case ES_REGULAR_DATA_AND_GOTO_NEXT_STEP:
		return ESCmd(u);
	}

	return ES_UNKNOWN;
}

static const char *ESCmd2str[] = {
	"0",
	"1",
	"REGULAR_DATA",				// 2
	"3",
	"GOTO_NEXT_STEP",			// 4
	"ABORT",				// 5
	"REGULAR_DATA_AND_GOTO_NEXT_STEP",	// 6
	0
};

const char *
ExtRemoteHandler::cmd2str(ESCmd cmd) const
{
	if (cmd == ES_UNKNOWN) {
		return "UNKNOWN";
	}

	return ESCmd2str[cmd];
}

ExtRemoteHandler::ExtRemoteHandler(mbsleep_t SleepTime, bool bReadForces, bool bLastReadForce)
: ExtFileHandlerBase(SleepTime, 0),
bReadForces(bReadForces), bLastReadForce(bLastReadForce)
{
	NO_OP;
}

ExtRemoteHandler::~ExtRemoteHandler(void)
{
    NO_OP;
}

void
ExtRemoteHandler::AfterPredict(void)
{
	bLastReadForce = false;
	bReadForces = true;
}

void
ExtRemoteHandler::Send_post(SendWhen when)
{
#if 0
	if (when == SEND_AFTER_CONVERGENCE) {
		bReadForces = true;
	}
#endif
	NO_OP;
}

bool
ExtRemoteHandler::ActOnCmd(uint8_t u)
{
	// TODO: act upon status code value
	switch (u2cmd(u)) {
	case ES_UNKNOWN:
		silent_cerr("ExtRemoteHandler: "
			"received unknown code (" << unsigned(u) << ")"
			<< std::endl);
		mbdyn_set_stop_at_end_of_iteration();
		return (bOK = false);

	case ES_REGULAR_DATA:
		break;

	case ES_GOTO_NEXT_STEP:
		// peer is done; do not read forces, keep using old
		bReadForces = false;
		return false;

	case ES_ABORT:
		silent_cout("ExtRemoteHandler: peer requested end of simulation"
			<< std::endl);
		mbdyn_set_stop_at_end_of_time_step();
		bReadForces = false;
		return false;

	case ES_REGULAR_DATA_AND_GOTO_NEXT_STEP:
		// peer is done; read forces for the last time
		bLastReadForce = true;
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return true;
}

bool
ExtRemoteHandler::Recv_post(void)
{
	if (bLastReadForce) {
		bReadForces = false;
	}

	return !bReadForces;
}

/* ExtSocketHandler moved to extsocket.h, extsocket.cc */

/* ExtFileHandlerEDGE moved to extedge.h, extedge.cc */

/* ExtForce - begin */

/*
 * Communication patterns:

	- loose coupling:
		- receive forces at first residual assembly
		- send motion after convergence

	- tight coupling:
		- send motion after predict
		- receive forces at each residual assembly
		- send motion at each update after solution

	- almost-tight coupling:
		- send motion after predict
		- receive forces every some residual assembly
		- send motion every some update after solution
		- send motion after convergence (to make sure
		  peer gets the last solution)
 */

/* Costruttore */
ExtForce::ExtForce(unsigned int uL,
	DataManager *pDM,
	ExtFileHandlerBase *pEFH,
	bool bSendAfterPredict,
	int iCoupling,
	flag fOut)
: Elem(uL, fOut), 
Force(uL, fOut),
c(iCoupling > COUPLING_LOOSE ? pDM : NULL),
pEFH(pEFH),
bSendAfterPredict(bSendAfterPredict),
iCoupling(iCoupling),
iCouplingCounter(0),
bFirstSend(true),
bFirstRecv(true)
{
	NO_OP;
}

ExtForce::~ExtForce(void)
{
	if (pEFH) {
		SAFEDELETE(pEFH);
	}
}

void
ExtForce::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints* h)
{
	bFirstSend = true;
	bFirstRecv = true;

	bool ok = false;
	if (pEFH->Prepare_pre()) {
		ok = Prepare(pEFH);
	}
	pEFH->Prepare_post(ok);
	if (!ok) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (iCoupling < COUPLING_TIGHT) {
		Send(ExtFileHandlerBase::SEND_AFTER_CONVERGENCE);
	}
}

void
ExtForce::Update(const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	/* If running tight coupling, send kinematics every iteration */
	/* NOTE: tight coupling may need relaxation */
	if ((iCoupling >= COUPLING_TIGHT && !((++iCouplingCounter)%iCoupling))
		|| (iCoupling == COUPLING_LOOSE && bSendAfterPredict))
	{
		ExtFileHandlerBase::SendWhen when;
		if (iCoupling == COUPLING_LOOSE) {
			when = ExtFileHandlerBase::SEND_AFTER_CONVERGENCE;
		} else if (bFirstSend) {
			when = ExtFileHandlerBase::SEND_FIRST_TIME;
		} else {
			when = ExtFileHandlerBase::SEND_REGULAR;
		}

		Send(when);
	}
}

/*
 * Elaborazione stato interno dopo la convergenza
 */
void
ExtForce::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	/* After prediction, mark next residual as first */
	iCouplingCounter = 0;
	bFirstSend = true;
	bFirstRecv = true;

	pEFH->AfterPredict();

	/* needed to send predicted data */
	if (bSendAfterPredict) {
		Update(X, XP);
	}
}

/*
 * Elaborazione stato interno dopo la convergenza
 */
void
ExtForce::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP)
{
	/* If not running tight coupling, send kinematics only at convergence */
	if (iCoupling != COUPLING_TIGHT && !bSendAfterPredict) {
		Send(ExtFileHandlerBase::SEND_AFTER_CONVERGENCE);
#if 0
		if (bRemoveIn) {
			Unlink();
		}
#endif
	}
}

/*
 * Send output to peer
 */
void
ExtForce::Send(ExtFileHandlerBase::SendWhen when)
{
	if (pEFH->Send_pre(when)) {
		Send(pEFH, when);
		bFirstSend = false;
	}
	pEFH->Send_post(when);
}

void
ExtForce::Recv(void)
{
	if ((iCoupling >= COUPLING_TIGHT && !bFirstSend && !(iCouplingCounter%iCoupling))
		|| ((iCoupling == COUPLING_LOOSE || iCoupling == COUPLING_STAGGERED) && bFirstRecv))
	{
		if (pEFH->Recv_pre()) {
			Recv(pEFH);
			bFirstRecv = false;
		}

		if (pEFH->Recv_post()) {
			// TODO: need to handle requests for end of simulation
			c.Set(Converged::CONVERGED);
		}
	}

	if (iCoupling == COUPLING_NONE) {
		c.Set(Converged::CONVERGED);
	}
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

void
ReadExtFileParams(DataManager* pDM,
	MBDynParser& HP, 
	unsigned int uLabel,
	mbsleep_t& SleepTime,
	std::streamsize& Precision)
{
	mbsleep_t MinSleepTime = SleepTime;
	if (HP.IsKeyWord("sleep" "time")) {
		SleepTime = HP.GetTimeout(SleepTime);
		if (SleepTime < MinSleepTime ) {
			silent_cerr("ExtForce(" << uLabel << "): "
				"invalid sleep time " << SleepTime
				<< " less than " << MinSleepTime
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	Precision = 6;
	if (HP.IsKeyWord("precision")) {
		if (!HP.IsKeyWord("default")) {
			Precision = HP.GetInt();
			if (Precision < 0) {
				silent_cerr("ExtForce(" << uLabel << "): "
					"invalid precision value "
					"\"" << Precision << "\""
					" at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}
}


static ExtFileHandlerBase *
ReadExtFileHandler(DataManager* pDM,
	MBDynParser& HP, 
	unsigned int uLabel)
{
	if (HP.IsKeyWord("EDGE")) {
		return ReadExtFileHandlerEDGE(pDM, HP, uLabel);

	} else if (HP.IsKeyWord("socket")) {
		return ReadExtSocketHandler(pDM, HP, uLabel);

	}
#ifdef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP
	else if (HP.IsKeyWord("shared" "memory")) {
		return ReadExtSharedMemHandler(pDM, HP, uLabel);

	}
#endif // HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP
	else {
	    silent_cerr("ExtForce(" << uLabel << "): "
			"unrecognised communicator type "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ExtFileHandlerBase *pEFH = 0;

	// default
	const char	*s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForce(" << uLabel << "): "
			"unable to get input file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	std::string fout = s;

	bool bNoClobberOut = false;
	if (HP.IsKeyWord("no" "clobber")) {
		bNoClobberOut = true;
	}

	mbsleep_t SleepTime = mbsleep_init(1);
	std::streamsize Precision = 0;
	ReadExtFileParams(pDM, HP, uLabel, SleepTime, Precision);

	SAFENEWWITHCONSTRUCTOR(pEFH, ExtFileHandler,
		ExtFileHandler(fin, bUnlinkIn, fout, bNoClobberOut,
			SleepTime, Precision));

	return pEFH;
}

void
ReadExtForce(DataManager* pDM, 
	MBDynParser& HP, 
	unsigned int uLabel,
	ExtFileHandlerBase*& pEFH,
	bool& bSendAfterPredict,
	int& iCoupling)
{
	pEFH = ReadExtFileHandler(pDM, HP, uLabel);

	iCoupling = ExtForce::COUPLING_LOOSE;
	if (HP.IsKeyWord("coupling")) {
		if (HP.IsKeyWord("none")) {
			iCoupling = ExtForce::COUPLING_NONE;

		} else if (HP.IsKeyWord("staggered")) {
			iCoupling = ExtForce::COUPLING_STAGGERED;

		} else if (HP.IsKeyWord("loose")) {
			iCoupling = ExtForce::COUPLING_LOOSE;

		} else if (HP.IsKeyWord("tight")) {
			iCoupling = ExtForce::COUPLING_TIGHT;

		} else {
			iCoupling = HP.GetInt();
			if (iCoupling < -1) {
				silent_cerr("ExtForce(" << uLabel << "): "
					"invalid coupling value "
					"\"" << iCoupling << "\""
					" at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	if (iCoupling >= ExtForce::COUPLING_LOOSE) {
		bSendAfterPredict = true;
		if (HP.IsKeyWord("send" "after" "predict")) {
			if (!HP.GetYesNo(bSendAfterPredict)) {
				silent_cerr("ExtForce(" << uLabel << "): "
					"invalud \"send after predict\" value "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

	} else {
		bSendAfterPredict = false;
	}
}
