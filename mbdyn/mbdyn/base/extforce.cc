/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2015
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
#include "except.h"
#include "solver.h"

#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <cstdlib>
#include <unistd.h>
#include <cerrno>

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
		int save_errno = errno;
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

/* ExtSocketHandler - begin */

#ifdef USE_SOCKET

ESCmd
ExtSocketHandler::u2cmd(unsigned u) const
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
ExtSocketHandler::cmd2str(ESCmd cmd) const
{
	if (cmd == ES_UNKNOWN) {
		return "UNKNOWN";
	}

	return ESCmd2str[cmd];
}

ExtSocketHandler::ExtSocketHandler(UseSocket *pUS, mbsleep_t SleepTime,
	int recv_flags, int send_flags)
: ExtFileHandlerBase(SleepTime, 0),
pUS(pUS), recv_flags(recv_flags), send_flags(send_flags),
bReadForces(true), bLastReadForce(false)
{
	NO_OP;
}

ExtSocketHandler::~ExtSocketHandler(void)
{
	if (pUS->GetSock() >= 0) {
		uint8_t u = ES_ABORT;

		// ignore result
		(void)send(pUS->GetSock(), (void *)&u, sizeof(u), send_flags);
	}
	SAFEDELETE(pUS);
}

ExtFileHandlerBase::Negotiate
ExtSocketHandler::NegotiateRequest(void) const
{
	if (pUS->Create()) {
		return ExtFileHandlerBase::NEGOTIATE_SERVER;

	} else {
		return ExtFileHandlerBase::NEGOTIATE_CLIENT;
	}
}

bool
ExtSocketHandler::Prepare_pre(void)
{
	uint8_t u;
	ssize_t rc;

	switch (NegotiateRequest()) {
	case ExtFileHandlerBase::NEGOTIATE_CLIENT:
		u = ES_NEGOTIATION;

		rc = send(pUS->GetSock(), (void *)&u, sizeof(u), send_flags);
		if (rc == -1) {
			int save_errno = errno;
			silent_cerr("ExtSocketHandler: negotiation request send() failed "
				"(" << save_errno << ": " << strerror(save_errno) << ")"
				<< std::endl);
			return (bOK = false);

		} else if (rc != sizeof(u)) {
			silent_cerr("ExtSocketHandler: negotiation request send() failed "
				"(sent " << rc << " bytes "
				"instead of " << sizeof(u) << ")"
				<< std::endl);
			return (bOK = false);
		}
		break;

	case ExtFileHandlerBase::NEGOTIATE_SERVER:
		rc = recv(pUS->GetSock(), (void *)&u, sizeof(u), recv_flags);
		if (rc == -1) {
			return (bOK = false);

		} else if (rc != sizeof(u)) {
			return (bOK = false);
		}

		if (u != ES_NEGOTIATION) {
			silent_cerr("ExtSocketHandler: negotiation request recv() failed"
				<< std::endl);
			return (bOK = false);
		}
		break;

	default:
		ASSERT(0);
	}

	return true;
}

void
ExtSocketHandler::Prepare_post(bool ok)
{
	if (NegotiateRequest()) {
		uint8_t u = ok ? ES_OK : ES_ABORT;

		ssize_t rc = send(pUS->GetSock(), (void *)&u, sizeof(u),
			send_flags);
		if (rc == -1) {
			int save_errno = errno;
			silent_cerr("ExtSocketHandler: negotiation response send() failed "
				"(" << save_errno << ": " << strerror(save_errno) << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (rc != sizeof(u)) {
			silent_cerr("ExtSocketHandler: negotiation response send() failed "
				"(sent " << rc << " bytes "
				"instead of " << sizeof(u) << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		uint8_t u;
		ssize_t rc = recv(pUS->GetSock(), (void *)&u, sizeof(u),
			recv_flags);
		if (rc == -1) {
			int save_errno = errno;
			silent_cerr("ExtSocketHandler: negotiation response recv() failed "
				"(" << save_errno << ": " << strerror(save_errno) << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (rc != sizeof(u)) {
			silent_cerr("ExtSocketHandler: negotiation response recv() failed "
				"(received " << rc << " bytes "
				"instead of " << sizeof(u) << ")"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		}

		if (u != ES_OK) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
}

void
ExtSocketHandler::AfterPredict(void)
{
	bLastReadForce = false;
	bReadForces = true;
}

bool
ExtSocketHandler::Send_pre(SendWhen when)
{
	if (!bReadForces || !bOK) {
		return false;
	}

	uint8_t u;
	if (when == SEND_AFTER_CONVERGENCE) {
		u = ES_REGULAR_DATA_AND_GOTO_NEXT_STEP;
	} else {
		u = ES_REGULAR_DATA;
	}
	ssize_t rc = send(pUS->GetSock(), (void *)&u, sizeof(u), send_flags);
	if (rc == -1) {
		int save_errno = errno;
		silent_cerr("ExtSocketHandler: send() failed "
			"(" << save_errno << ": " << strerror(save_errno) << ")"
			<< std::endl);
		mbdyn_set_stop_at_end_of_iteration();
		return (bOK = false);

	} else if (rc != sizeof(u)) {
		silent_cerr("ExtSocketHandler: send() failed "
			"(sent " << rc << " bytes "
			"instead of " << sizeof(u) << ")"
			<< std::endl);
		mbdyn_set_stop_at_end_of_iteration();
		return (bOK = false);
	}

	return true;
}

void
ExtSocketHandler::Send_post(SendWhen when)
{
#if 0
	if (when == SEND_AFTER_CONVERGENCE) {
		bReadForces = true;
	}
#endif
	NO_OP;
}

bool
ExtSocketHandler::Recv_pre(void)
{
	if (!bReadForces) {
		return false;
	}

	uint8_t u = 0;

	if (SleepTime != 0) {
		for ( ; ; ) {
			ssize_t rc;

			rc = recv(pUS->GetSock(), (void *)&u, sizeof(u),
				recv_flags | MSG_DONTWAIT);
			if (rc != -1) {
				break;
			}

			int save_errno = errno;
	
			if (errno != EAGAIN) {
				silent_cerr("ExtSocketHandler: "
					"recv() failed (" << save_errno << ": "
					<< strerror(save_errno) << ")"
					<< std::endl);
				mbdyn_set_stop_at_end_of_iteration();
				return (bOK = false);
			}

			if (mbdyn_stop_at_end_of_iteration()) {
				return (bOK = false);
			}

			mbsleep(&SleepTime);
		}

	} else {
		// wait until the status code is returned
		ssize_t rc = recv(pUS->GetSock(), (void *)&u, sizeof(u), recv_flags);
		if (rc == -1) {
			int save_errno = errno;
	
			if (errno != EAGAIN) {
				silent_cerr("ExtSocketHandler: "
					"recv() failed (" << save_errno << ": "
					<< strerror(save_errno) << ")"
					<< std::endl);
				mbdyn_set_stop_at_end_of_iteration();
				return (bOK = false);
			}

		} else if (rc != sizeof(u)) {
			silent_cerr("ExtSocketHandler: "
				"recv()=" << rc << " (expected " << sizeof(u) << ")" 
				<< std::endl);
			mbdyn_set_stop_at_end_of_iteration();
			return (bOK = false);
		}
	}

	// TODO: act upon status code value
	switch (u2cmd(u)) {
	case ES_UNKNOWN:
		silent_cerr("ExtSocketHandler: "
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
		silent_cout("ExtSocketHandler: peer requested end of simulation"
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
ExtSocketHandler::Recv_post(void)
{
	if (bLastReadForce) {
		bReadForces = false;
	}

	return !bReadForces;
}

int
ExtSocketHandler::GetOutFileDes(void)
{
	return pUS->GetSock();
}

int
ExtSocketHandler::GetSendFlags(void) const
{
	return send_flags;
}

int
ExtSocketHandler::GetInFileDes(void)
{
	return pUS->GetSock();
}

int
ExtSocketHandler::GetRecvFlags(void) const
{
	return recv_flags;
}

#endif // USE_SOCKET

/* ExtSocketHandler - end */

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
ReadExtSocketHandler(DataManager* pDM,
	MBDynParser& HP, 
	unsigned int uLabel)
{
#ifdef USE_SOCKET
	ExtFileHandlerBase *pEFH = 0;

	bool bCreate = false;
	unsigned short int port = (unsigned short int)-1; 
	std::string host;
	std::string path;

	if (HP.IsKeyWord("create")) {
		if (!HP.GetYesNo(bCreate)) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"\"create\" must be either \"yes\" or \"no\" "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
		
	if (HP.IsKeyWord("local") || HP.IsKeyWord("path")) {
		const char *m = HP.GetFileName();
		
		if (m == 0) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"unable to read local path "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		
		path = m;
	}
	
	if (HP.IsKeyWord("port")) {
		if (!path.empty()) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"cannot specify port "
				"for a local socket "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);		
		}

		int p = HP.GetInt();
		/* Da sistemare da qui */
#ifdef IPPORT_USERRESERVED
		if (p <= IPPORT_USERRESERVED) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"cannot listen on reserved port "
				<< port << ": less than "
				"IPPORT_USERRESERVED=" << IPPORT_USERRESERVED
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		/* if #undef'd, don't bother checking;
		 * the OS will do it for us */
#endif /* IPPORT_USERRESERVED */

		port = p;
	}

	if (HP.IsKeyWord("host")) {
		if (!path.empty()) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"cannot specify host for a local socket "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);		
		}

		const char *h;
		
		h = HP.GetStringWithDelims();
		if (h == 0) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"unable to read host "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		host = h;

	} else if (path.empty() && !bCreate) {
		silent_cerr("ExtSocketHandler"
			"(" << uLabel << "): "
			"host undefined "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	int socket_type = SOCK_STREAM;
	if (HP.IsKeyWord("socket" "type")) {
		if (HP.IsKeyWord("udp")) {
			socket_type = SOCK_DGRAM;

		} else if (!HP.IsKeyWord("tcp")) {
			silent_cerr("ExtSocketHandler(" << uLabel << "\"): "
				"invalid socket type "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	if ((socket_type == SOCK_DGRAM) && !bCreate) {
		silent_cerr("ExtSocketHandler(" << uLabel << "\"): "
			"socket type=upd incompatible with create=no "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// we want to block until the whole chunk is received
	int recv_flags = 0;
	int send_flags = 0;
#ifdef MSG_WAITALL
	recv_flags |= MSG_WAITALL;
#endif // MSG_WAITALL

	while (HP.IsArg()) {
		if (HP.IsKeyWord("signal")) {
#ifdef MSG_NOSIGNAL
			recv_flags &= ~MSG_NOSIGNAL;
			send_flags &= ~MSG_NOSIGNAL;
#else // ! MSG_NOSIGNAL
			silent_cout("ExtSocketHandler"
				"(" << uLabel << "): "
				"MSG_NOSIGNAL not defined (ignored) "
				"at line " << HP.GetLineData()
				<< std::endl);
#endif // ! MSG_NOSIGNAL

		// not honored by recv(2)
		} else if (HP.IsKeyWord("no" "signal")) {
#ifdef MSG_NOSIGNAL
			recv_flags |= MSG_NOSIGNAL;
			send_flags |= MSG_NOSIGNAL;
#else // ! MSG_NOSIGNAL
			silent_cout("ExtSocketHandler"
				"(" << uLabel << "): "
				"MSG_NOSIGNAL not defined (ignored) "
				"at line " << HP.GetLineData()
				<< std::endl);
#endif // ! MSG_NOSIGNAL

		} else {
			break;
		}
	}

	UseSocket *pUS = 0;
	if (path.empty()) {
		if (port == (unsigned short int)(-1)) {
			silent_cerr("ExtSocketHandler"
				"(" << uLabel << "): "
				"port missing"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		SAFENEWWITHCONSTRUCTOR(pUS, UseInetSocket, UseInetSocket(host, port, socket_type, bCreate));

	} else {
		SAFENEWWITHCONSTRUCTOR(pUS, UseLocalSocket, UseLocalSocket(path, socket_type, bCreate));
	}

	if (bCreate) {
		pDM->RegisterSocketUser(pUS);

	} else {
		pUS->Connect();
	}

	mbsleep_t SleepTime = mbsleep_init(0);
	std::streamsize Precision = 0;
	ReadExtFileParams(pDM, HP, uLabel, SleepTime, Precision);
	// NOTE: so far, precision is ignored

	SAFENEWWITHCONSTRUCTOR(pEFH, ExtSocketHandler,
		ExtSocketHandler(pUS, SleepTime, recv_flags, send_flags));

	return pEFH;
#else // ! USE_SOCKET
	silent_cerr("ExtSocketHandler not supported" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_SOCKET
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
