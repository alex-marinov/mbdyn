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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cerrno>
#include "myassert.h"
#include "solver.h"
#include "solver_impl.h"
#include "rtsolver.h"
#include "rtposixsolver.h"
 
#ifdef USE_RT

/* RTPOSIXSolver - begin */

RTPOSIXSolver::RTPOSIXSolver(Solver *pS,
	RTMode eRTMode,
	unsigned long lRTPeriod,
	unsigned long RTStackSize,
	bool bRTAllowNonRoot,
	int RTCpuMap,
	bool bNoOutput)
: RTSolverBase(pS, eRTMode, lRTPeriod, RTStackSize, bRTAllowNonRoot, RTCpuMap, bNoOutput),
clock_flags(TIMER_ABSTIME)
{
	NO_OP;
}

RTPOSIXSolver::~RTPOSIXSolver(void)
{
	NO_OP;
}

// write contribution to restart file
std::ostream&
RTPOSIXSolver::Restart(std::ostream& out) const
{
	return out;
}

// very first setup, to be always performed
void
RTPOSIXSolver::Setup(void)
{
	// allow-nonroot should go here; something like
	// capset(CAP_SYS_NICE);

	struct sched_param sched;
	int policy = SCHED_FIFO;
	int priority = 1;
	int min_priority = sched_get_priority_min(policy);

	sched.sched_priority = sched_get_priority_max(policy) - priority;

	if (sched.sched_priority < min_priority) {
		sched.sched_priority = min_priority;
	}

        if (sched_setscheduler(0, policy, &sched) < 0) {
		silent_cerr("RTPOSIXSolver: sched_setscheduler failed" << std::endl);
		if (bRTAllowNonRoot) {
			// FIXME: hack (let it run without priority)
			return;
		}
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

#ifdef HAVE_SCHED_SETAFFINITY
	if (RTCpuMap != 0xFF) {
		cpu_set_t cpuset;

		CPU_ZERO(&cpuset);
		for (int cpu = 0; cpu < 8; cpu++) {
			if ((RTCpuMap >> cpu) & 0x1) {
				CPU_SET(cpu, &cpuset);
			}
		}

		if (sched_setaffinity(0, sizeof(cpu_set_t), &cpuset)) {
			silent_cerr("RTPOSIXSolver: sched_setaffinity failed" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
#endif // HAVE_SCHED_SETAFFINITY
}

// to be performed when stop is commanded by someone else
void
RTPOSIXSolver::StopCommanded(void)
{
	// do nothing
	NO_OP;
}

// write real-time related message when stop commanded by someone else
void
RTPOSIXSolver::Log(void)
{
	// TODO: print info about overruns?
	NO_OP;
}

// wait for period to expire
void
RTPOSIXSolver::Wait(void)
{
	if (RTWaitPeriod()) {
		if (RTSteps == 0) {
			clock_gettime(CLOCK_MONOTONIC, &t0);
			t = t0;
		}

		t.tv_nsec += lRTPeriod;
		if (t.tv_nsec >= 1000000000L) {
			t.tv_nsec -= 1000000000L;
			t.tv_sec++;
		}

		int rc = clock_nanosleep(CLOCK_MONOTONIC, clock_flags, &t, NULL);
		switch (rc) {
		case 0:
			break;

		case EINTR:
			silent_cerr("RTPOSIXSolver: clock_nanosleep interrupted" << std::endl);
			break;

		default:
			silent_cerr("RTPOSIXSolver: clock_nanosleep failed "
				"(rc=" << rc << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

#if 0
	} else if (RTSemWait()) {
#endif
	} /* else RTBlockingIO(): do nothing */

	RTSteps++;
}

#endif // USE_RT

RTSolverBase *
ReadRTPOSIXSolver(Solver *pS, MBDynParser& HP)
{
#ifdef USE_RT
	RTSolverBase::RTMode eRTMode;
	unsigned long lRTPeriod;
	unsigned long RTStackSize;
	bool bRTAllowNonRoot;
	int RTCpuMap;
	ReadRTParams(pS, HP, eRTMode, lRTPeriod, RTStackSize, bRTAllowNonRoot, RTCpuMap);

	bool bNoOutput = true;
	if (HP.IsKeyWord("output")) {
		bNoOutput = HP.GetYesNoOrBool(bNoOutput);
	}

	RTSolverBase *pRTSolver(0);
	SAFENEWWITHCONSTRUCTOR(pRTSolver, RTPOSIXSolver,
		RTPOSIXSolver(pS, eRTMode, lRTPeriod,
			RTStackSize, bRTAllowNonRoot, RTCpuMap,
			bNoOutput));
	return pRTSolver;

#else // !USE_RT
	silent_cerr("ReadRTPOSIXSolver: need to configure --with-rt "
		"to use POSIX realtime" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	return 0;
#endif // USE_RT
}

/* RTPOSIXSolver - end */

