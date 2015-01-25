/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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


#include "mbdefs.h"
#include "solver.h"
#include "solver_impl.h"
#include "rtaisolver.h"

#include "mbrtai_utils.h"

// RTAI log message
struct mbrtai_msg_t {
	int step;
	int time;
} msg;

/* RTAISolver - begin */

RTAISolver::RTAISolver(Solver *pS,
	RTMode eRTMode,
	long long lRTPeriod,
	unsigned long RTStackSize,
	bool bRTAllowNonRoot,
	int RTCpuMap,
	bool bRTHard,
	bool bRTlog,
	const std::string& LogProcName)
: RTSolverBase(pS, eRTMode, lRTPeriod,
	RTStackSize, bRTAllowNonRoot, RTCpuMap),
bRTHard(bRTHard),
bRTlog(bRTlog),
LogProcName(LogProcName),
lRTPeriod(lRTPeriod),
RTSemPtr_in(0),
RTSemPtr_out(0),
mbxlog(0),
RTStpFlag(0),
t_tot(0),
t0(0),
t1(0),
or_counter(0)
{
	ASSERT(lRTPeriod > 0);
	ASSERT(!bRTlog || !LogProcName.empty());
}

RTAISolver::~RTAISolver(void)
{
	if (mbxlog) {
		rtmbdyn_rt_mbx_delete(&mbxlog);
		mbxlog = 0;
	}
}

// write contribution to restart file
std::ostream&
RTAISolver::Restart(std::ostream& out) const
{
	out << "RTAI";

	out << ", reserve stack, " << RTStackSize;

	out << ", mode, ";
	switch (eRTMode) {
	case RTAISolver::MBRT_WAITPERIOD:
		out << "period, time step, " << lRTPeriod;
		break;

	case RTAISolver::MBRT_SEMAPHORE:
		out << "semaphore";
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (bRTAllowNonRoot) {
		out << ", allow nonroot";
	}

	if (RTCpuMap != 0xFF) {
		out << ", cpu map, " << RTCpuMap;
	}

	if (bRTHard) {
		out << ", hard realtime";
	}

	if (bRTlog) {
		out << ", realtime log, \"" << LogProcName << "\"";
	}

	return out;
}

// very first setup, to be always performed
void
RTAISolver::Setup(void)
{
	if (bRTAllowNonRoot) {
		rtmbdyn_rt_allow_nonroot_hrt();
	}

	ASSERT(::rtmbdyn_rtai_task == 0);

	/* Init RTAI; if init'ed, it will be shut down at exit */
	if (rtmbdyn_rt_task_init("MBDTSK", 1, 0, 0, RTCpuMap,
		&::rtmbdyn_rtai_task))
	{
		silent_cerr("RTAISolver: unable to init RTAI task" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

// initialization to be performed only if real-time is requested
void
RTAISolver::Init(void)
{
	/* Need timer */
	if (!rtmbdyn_rt_is_hard_timer_running()) {
		/* FIXME: ??? */
		silent_cout("RTAISolver: Hard timer is started by MBDyn"
			<< std::endl);
		rtmbdyn_rt_set_oneshot_mode();
		rtmbdyn_start_rt_timer(rtmbdyn_nano2count(1000000));
	}

	/*
	 * MBDyn can work in two ways:
	 * - internal timer
	 * - scheduled by an external signal
	 * only the first case is currently implemented
	 */
	if (RTWaitPeriod()) {
		long long t = rtmbdyn_rt_get_time();
		int r;

		/* Timer should be init'ed */
		ASSERT(t > 0);

		silent_cout("RTAISolver: Task: " << ::rtmbdyn_rtai_task
			<< "; time: " << t
			<< "; period: " << lRTPeriod
			<< std::endl);

		// NOTE: the period was in nanoseconds until now.
		lRTPeriod = rtmbdyn_nano2count(lRTPeriod);

		r = rtmbdyn_rt_task_make_periodic(::rtmbdyn_rtai_task,
			t, lRTPeriod);

		if (r) {
			silent_cerr("RTAISolver: "
				"rtmbdyn_rt_task_make_periodic() failed "
				"(" << r << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

#if 0
	} else if (RTSemWait()) {
		int r;

		/* FIXME: check args
		 * name should be configurable?
		 * initial value 0: non-blocking
		 */
		r = rtmbdyn_rt_sem_init("MBDSMI", 0, &RTSemPtr_in);
		if (r) {
			silent_cerr("rt_sem_init() failed ("
				<< r << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/* FIXME: check args
		 * name should be configurable?
		 * initial value 0: non-blocking
		 */
		r = rtmbdyn_rt_sem_init("MBDSMO", 0, &RTSemPtr_out);
		if (r) {
			silent_cerr("rt_sem_init() failed ("
				<< r << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
#endif
	}

	/* FIXME: should check whether RTStackSize is correctly set? */
	if (bRTlog) {
		silent_cout("RTAISolver: MBDyn starts overruns monitor "
			"(proc: \"" << LogProcName << "\")"
			<< std::endl);

		const char *mbxlogname = "logmb";
		if (rtmbdyn_rt_mbx_init(mbxlogname, sizeof(msg)*16, &mbxlog)) {
			bRTlog = false;
			silent_cerr("RTAISolver: cannot init log mailbox "
				"\"" << mbxlogname << "\""
				<< std::endl);

		} else {
			const char *nonroot =
				bRTAllowNonRoot ? "TRUE" : "FALSE";

			switch (fork()) {
			case 0: {
				char LogCpuMap[] = "0xFF";

				if (RTCpuMap != 0xFF) {
					/* MBDyn can use any cpu
					 * The overruns monitor will use any free cpu */
					snprintf(LogCpuMap, sizeof(LogCpuMap),
						"0x%02X", ~RTCpuMap);
				}

				if (strcmp(LogProcName.c_str(), "logproc") != 0) {
					if (execl(LogProcName.c_str(), LogProcName.c_str(),
						"MBDTSK", mbxlogname,
						LogCpuMap, nonroot, NULL) == 0)
					{
						break;
					}

					/* error */
					silent_cout("RTAISolver: cannot start "
						"log procedure "
						"\"" << LogProcName << "\"; "
						"using default" << std::endl);
				}

#ifdef HAVE_SETENV
				/* sets new path */
				/* BINPATH is the ${bindir} variable
				 * at configure time, defined in
				 * include/mbdefs.h.in */
				char *origpath = getenv("PATH");
				if (origpath == 0) {
					/* ?!? */
					setenv("PATH", ".:" BINPATH, 1);

				} else {
					std::string newpath = ".:" BINPATH ":";
					newpath += origpath;
					setenv("PATH", newpath.c_str(), 1);
				}
#endif // HAVE_SETENV

				/* start logger */
				if (execlp("logproc", "logproc", "MBDTSK",
			               	mbxlogname, LogCpuMap, nonroot, NULL)
					== -1)
				{
					silent_cout("RTAISolver: cannot start default "
						"log procedure \"logproc\""
						<< std::endl);
					/* FIXME: better give up logging? */
					bRTlog = false;
				}
				break;
			}

			case -1:
				silent_cerr("Cannot init log procedure" << std::endl);
				bRTlog = false;
				break;

			default:
				rtmbdyn_rt_sleep(rtmbdyn_nano2count(1000000000));
				break;
			}
		}
	}

	RTSolverBase::Init();
}

// check whether stop is commanded by real-time
bool
RTAISolver::IsStopCommanded(void)
{
	if (RTStpFlag == 1) {
		StopCommanded();
	}

	return (RTStpFlag != 0);
}

// to be performed when stop is commanded by someone else
void
RTAISolver::StopCommanded(void)
{
	if (bRTHard) {
		rtmbdyn_rt_make_soft_real_time();
	}
}

// write real-time related message when stop commanded by someone else
void
RTAISolver::Log(void)
{
	silent_cout("total overruns: " << or_counter  << std::endl
		  << "total overrun time: " << t_tot << " micros" << std::endl);
}

// wait for period to expire
void
RTAISolver::Wait(void)
{
	rtmbdyn_rt_receive_if(NULL, &RTStpFlag);

	t1 = rtmbdyn_rt_get_time();
	if (RTSteps >= 2 && t1 > (t0 + lRTPeriod)) {
		or_counter++;
		t_tot = t_tot + rtmbdyn_count2nano(t1 - t0 - lRTPeriod)/1000;

		if (bRTlog) {
			msg.step = RTSteps;
			msg.time = (int)rtmbdyn_count2nano(t1 - t0 - lRTPeriod)/1000;

			rtmbdyn_RT_mbx_send_if(0, 0, mbxlog, &msg, sizeof(msg));
		}
	}

	if (RTWaitPeriod()) {
		rtmbdyn_rt_task_wait_period();

#if 0
	} else if (RTSemWait()) {
		/* FIXME: semaphore must be configurable */
		mbdyn_rt_sem_wait(RTSemPtr_in);
#endif
	} /* else RTBlockingIO(): do nothing */

	t0 = rtmbdyn_rt_get_time();

	if (RTSteps == 2 && bRTHard) {
		/* make hard real time */
		rtmbdyn_rt_make_hard_real_time();
	}

	RTSteps++;
}

/* RTAISolver - end */

RTSolverBase *
ReadRTAISolver(Solver *pS, MBDynParser& HP)
{
	RTSolverBase::RTMode eRTMode;
	unsigned long lRTPeriod;
	unsigned long RTStackSize;
	bool bRTAllowNonRoot;
	int RTCpuMap;
	ReadRTParams(pS, HP, eRTMode, lRTPeriod, RTStackSize, bRTAllowNonRoot, RTCpuMap);

	bool bRTHard = false;
	if (HP.IsKeyWord("hard" "real" "time")) {
		bRTHard = true;
	}

	bool bRTlog(false);
	std::string LogProcName;
	if (HP.IsKeyWord("real" "time" "log")) {
		if (HP.IsKeyWord("file" "name")){
			const char *m = HP.GetFileName();
			if (m == 0) {
				silent_cerr("RTAISolver: unable to get "
					"log process name (\"file name\") "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			LogProcName = m;

		} else {
			// built-in log process
			LogProcName = "logproc";
		}

		bRTlog = true;
	}

	RTSolverBase *pRTSolver(0);
	SAFENEWWITHCONSTRUCTOR(pRTSolver, RTAISolver,
		RTAISolver(pS, eRTMode, lRTPeriod,
			RTStackSize, bRTAllowNonRoot, RTCpuMap,
			bRTHard, bRTlog, LogProcName));

	return pRTSolver;
}

