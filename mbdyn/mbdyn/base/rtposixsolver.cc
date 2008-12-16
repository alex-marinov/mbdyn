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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

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
	int RTCpuMap)
: RTSolverBase(pS, eRTMode, lRTPeriod, RTStackSize, bRTAllowNonRoot, RTCpuMap),
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
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
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
	if (RTSteps == 0) {
		clock_gettime(CLOCK_MONOTONIC, &t0);
		t = t0;
	}

	t.tv_nsec += lRTPeriod;
	if (t.tv_nsec >= 1000000000) {
		t.tv_nsec -= 1000000000;
		t.tv_sec++;
	}

	switch (clock_nanosleep(CLOCK_MONOTONIC, clock_flags, &t, NULL)) {
	case 0:
		break;

	case -EINVAL:
	case -EFAULT:
		silent_cerr("RTSolver: clock_nanosleep failed" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

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

	RTSolverBase *pRTSolver(0);
	SAFENEWWITHCONSTRUCTOR(pRTSolver, RTPOSIXSolver,
		RTPOSIXSolver(pS, eRTMode, lRTPeriod,
			RTStackSize, bRTAllowNonRoot, RTCpuMap));
	return pRTSolver;

#else // !USE_RT
	silent_cerr("RTPOSIXSolver: need to configure --with-rt "
		"to use POSIX realtime" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	return 0;
#endif // USE_RT
}

/* RTPOSIXSolver - end */

