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
 
#ifndef RTAISOLVER_H
#define RTAISOLVER_H  

#include "rtsolver.h"

#ifdef USE_RTAI

/* RTAISolver - begin */

class RTAISolver : public RTSolverBase {
public:
	enum RTMode {
		MBRTAI_UNKNOWN,
		MBRTAI_WAITPERIOD,
		MBRTAI_SEMAPHORE,
		MBRTAI_LASTMODE
	};

private:
	bool bRTAllowNonRoot;
	RTMode eRTMode;
	bool bRTHard;

	/* if eRTMode == MBRTAI_WAITPERIOD */
	long long lRTPeriod;

	/* if eRTMode == MBRTAI_SEMAPHORE */
	void *RTSemPtr_in;
	void *RTSemPtr_out;

	int RTCpuMap;
	bool bRTlog;
	char *LogProcName;

	bool RTWaitPeriod(void) const {
		return (eRTMode == MBRTAI_WAITPERIOD);
	};

	bool RTSemaphore(void) const {
		return (eRTMode == MBRTAI_SEMAPHORE);
	};

	void *mbxlog;
	int RTStpFlag;
	volatile int RTSteps;
        int t_tot;
	long long t0, t1;
	int or_counter;

public:
	RTAISolver(Solver *pS,
		unsigned long RTStackSize,
		bool bRTAllowNonRoot,
		RTMode eRTMode,
		bool bRTHard,
		long long lRTPeriod,
		int RTCpuMap,
		bool bRTlog,
		char *LogProcName);
	~RTAISolver(void);

	// very first setup, to be always performed
	void Setup(void);
	// initialization to be performed only if real-time is requested
	void Init(void);
	// check whether stop is commanded by real-time
	bool IsStopCommanded(void);
	// to be performed when stop is commanded by someone else
	void StopCommanded(void);
	// write real-time related message when stop commanded by someone else
	void Log(void);
	// wait for period to expire
	void Wait(void);
};

/* RTAISolver - end */

#endif /* USE_RTAI */

extern RTSolverBase *
ReadRTAISolver(Solver *pS, MBDynParser& HP);

#endif // RTAISOLVER_H

