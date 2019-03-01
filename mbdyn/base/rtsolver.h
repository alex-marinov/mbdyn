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
 
#ifndef RTSOLVER_H
#define RTSOLVER_H  

/* RTSolverBase - begin */

class Solver;

class RTSolverBase {
public:
	enum RTMode {
		MBRT_UNKNOWN = -1,

		MBRT_WAITPERIOD,
		MBRT_SEMAPHORE,
		MBRT_IO,

		MBRT_LASTMODE
	};

protected:
	Solver *pS;

	RTMode eRTMode;

	/* if eRTMode == MBRT_WAITPERIOD */
	unsigned long lRTPeriod;

	unsigned long RTStackSize;
	bool bRTAllowNonRoot;
	int RTCpuMap;

	bool bNoOutput;

	bool RTWaitPeriod(void) const {
		return (eRTMode == MBRT_WAITPERIOD);
	};

	bool RTSemWait(void) const {
		return (eRTMode == MBRT_SEMAPHORE);
	};

	bool RTBlockingIO(void) const {
		return (eRTMode == MBRT_IO);
	};

	volatile int RTSteps;

public:
	RTSolverBase(Solver *pS,
		RTMode eRTMode,
		unsigned long lRTPeriod,
		unsigned long RTStackSize,
		bool bRTAllowNonRoot,
		int RTCpuMap,
		bool bNoOutput = true);
	virtual ~RTSolverBase(void);

	// write contribution to restart file
	virtual std::ostream& Restart(std::ostream& out) const = 0;
	// very first setup, to be always performed
	virtual void Setup(void) = 0;
	// initialization to be performed only if real-time is requested
	virtual void Init(void);
	// check whether stop is commanded by real-time
	virtual bool IsStopCommanded(void);
	// to be performed when stop is commanded by someone else
	virtual void StopCommanded(void) = 0;
	// write real-time related message when stop commanded by someone else
	virtual void Log(void) = 0;
	// wait for period to expire
	virtual void Wait(void) = 0;
};

/* RTSolverBase - end */

extern void
ReadRTParams(Solver *pS, MBDynParser& HP,
	RTSolverBase::RTMode& eRTMode,
	unsigned long& lRTPeriod,
	unsigned long& RTStackSize,
	bool& bRTAllowNonRoot,
	int& RTCpuMap);

extern RTSolverBase *
ReadRTSolver(Solver *pS, MBDynParser& HP);

#endif // RTSOLVER_H

