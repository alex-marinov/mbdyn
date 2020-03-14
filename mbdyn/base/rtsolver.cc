/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "myassert.h"
#include "solver.h"
#include "solver_impl.h"
#include "rtsolver.h"
#include "rtaisolver.h"
#include "rtposixsolver.h"
#include "ac/sys_sysinfo.h"
 
/* RTSolverBase - begin */

RTSolverBase::RTSolverBase(Solver *pS,
	RTMode eRTMode,
	unsigned long lRTPeriod,
	unsigned long RTStackSize,
	bool bRTAllowNonRoot,
	int RTCpuMap,
	bool bNoOutput)
: pS(pS),
eRTMode(eRTMode),
lRTPeriod(lRTPeriod),
RTStackSize(RTStackSize),
bRTAllowNonRoot(bRTAllowNonRoot),
RTCpuMap(RTCpuMap),
bNoOutput(bNoOutput)
{
	ASSERT(RTStackSize > 0);
	ASSERT(lRTPeriod > 0);
}

RTSolverBase::~RTSolverBase(void)
{
	NO_OP;
}

void
RTSolverBase::Init(void)
{
	/* if using real-time, clear out any type of output */
	if (bNoOutput) {
		pS->SetNoOutput();
	}
	RTSteps = 0;
	mbdyn_reserve_stack(RTStackSize);
}

// check whether stop is commanded by real-time
bool
RTSolverBase::IsStopCommanded(void)
{
	return false;
}

/* RTSolverBase - end */

void
ReadRTParams(Solver *pS, MBDynParser& HP,
	RTSolverBase::RTMode& eRTMode,
	unsigned long& lRTPeriod,
	unsigned long& RTStackSize,
	bool& bRTAllowNonRoot,
	int& RTCpuMap)
{
	eRTMode = RTSolverBase::MBRT_UNKNOWN;
	/* FIXME: use a safe default? */
	if (HP.IsKeyWord("mode")) {
		if (HP.IsKeyWord("period")) {
			eRTMode = RTSolverBase::MBRT_WAITPERIOD;

		} else if (HP.IsKeyWord("semaphore")) {
			/* FIXME: not implemented yet ... */
			eRTMode = RTSolverBase::MBRT_SEMAPHORE;

		} else if (HP.IsKeyWord("io")) {
			eRTMode = RTSolverBase::MBRT_IO;

		} else {
			silent_cerr("RTSolver: unknown realtime mode "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	switch (eRTMode) {
	case RTSolverBase::MBRT_UNKNOWN:
		eRTMode = RTSolverBase::MBRT_WAITPERIOD;
		silent_cerr("RTSolver: unknown realtime mode; "
			"assuming periodic "
			"at line " << HP.GetLineData()
			<< std::endl);
		// fallthru

	case RTSolverBase::MBRT_WAITPERIOD:
		lRTPeriod = -1;
		if (HP.IsKeyWord("time" "step")) {
			long long p = HP.GetInt();

			if (p <= 0) {
				silent_cerr("RTSolver: illegal time step "
					<< p << " at line "
					<< HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			lRTPeriod = p;

		} else {
			silent_cerr("RTSolver: need a time step for real time "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		break;

	case RTSolverBase::MBRT_SEMAPHORE:
		// impossible, right now
		break;

	case RTSolverBase::MBRT_IO:
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	RTStackSize = 1024;
	if (HP.IsKeyWord("reserve" "stack")) {
		long size = HP.GetInt();

		if (size <= 0) {
			silent_cerr("RTSolver: illegal stack size "
				<< size << " at line "
				<< HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		RTStackSize = size;
	}

	bRTAllowNonRoot = false;
	if (HP.IsKeyWord("allow" "nonroot")) {
		bRTAllowNonRoot = true;
	}

	RTCpuMap = 0xFF;
	if (HP.IsKeyWord("cpu" "map")) {
#ifndef HAVE_SCHED_SETAFFINITY
		silent_cerr("warning: \"cpu map\" unsupported (ignored)" << std::endl);
#endif // ! HAVE_SCHED_SETAFFINITY

		int cpumap = HP.GetInt();
		// NOTE: there is a hard limit at 4 CPU
		int ncpu = std::min(get_nprocs(), 4);
		int newcpumap = (2 << (ncpu - 1)) - 1;

		/* i bit non legati ad alcuna cpu sono posti a zero */
		newcpumap &= cpumap;
		if (newcpumap < 1 || newcpumap > 0xFF) {
			char buf[16];
			snprintf(buf, sizeof(buf), "0x%02X", (unsigned)cpumap);
			silent_cerr("RTSolver: illegal cpu map " << buf
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		RTCpuMap = newcpumap;
	}
}

RTSolverBase *
ReadRTSolver(Solver *pS, MBDynParser& HP)
{
	if (HP.IsKeyWord("POSIX")) {
		return ReadRTPOSIXSolver(pS, HP);
	}

#ifdef USE_RTAI
	if (!HP.IsKeyWord("RTAI")) {
		pedantic_cerr("ReadRTSolver: missing real-time model; "
			"assuming RTAI at line " << HP.GetLineData()
			<< std::endl);
	}

	return ReadRTAISolver(pS, HP);

#else // ! USE_RTAI
	if (HP.IsKeyWord("RTAI")) {
		silent_cerr("ReadRTSolver: need to configure --with-rtai "
			"to use RTAI real-time solver "
			"at line " << HP.GetLineData() << std::endl);

	} else {
		silent_cerr("ReadRTSolver: need to configure --with-rtai "
			"to use default RTAI real-time solver "
			"at line " << HP.GetLineData() << std::endl);
	}

	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif // ! USE_RTAI
}

/* RTSolver - end */

