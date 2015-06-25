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
 
#ifndef RTPOSIXSOLVER_H
#define RTPOSIXSOLVER_H  

#ifdef USE_RT

#include "rtsolver.h"

/* RTPOSIXSolver - begin */

class RTPOSIXSolver : public RTSolverBase {
private:
	int clock_flags;
	struct timespec t0, t;

public:
	RTPOSIXSolver(Solver *pS,
		RTMode eRTMode,
		unsigned long lRTPeriod,
		unsigned long RTStackSize,
		bool bRTAllowNonRoot,
		int RTCpuMap,
		bool bNoOutput);
	~RTPOSIXSolver(void);

	// write contribution to restart file
	std::ostream& Restart(std::ostream& out) const;
	// very first setup, to be always performed
	void Setup(void);
	// to be performed when stop is commanded by someone else
	void StopCommanded(void);
	// write real-time related message when stop commanded by someone else
	void Log(void);
	// wait for period to expire
	void Wait(void);
};

/* RTPOSIXSolver - end */

#endif // USE_RT

extern RTSolverBase *
ReadRTPOSIXSolver(Solver *pS, MBDynParser& HP);

#endif // RTPOSIXSOLVER_H

