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
#include "rtaisolver.h"
 
/* RTSolverBase - begin */

RTSolverBase::RTSolverBase(Solver *pS, unsigned long RTStackSize)
: pS(pS),
RTStackSize(RTStackSize)
{
	NO_OP;
}

RTSolverBase::~RTSolverBase(void)
{
	NO_OP;
}

void
RTSolverBase::Init(void)
{
	/* if using real-time, clear out any type of output */
	pS->SetNoOutput();

	ASSERT(RTStackSize > 0);
	mbdyn_reserve_stack(RTStackSize);
}

// check whether stop is commanded by real-time
bool
RTSolverBase::IsStopCommanded(void)
{
	return false;
}

/* RTSolverBase - end */

/* RTSolver - begin */

RTSolver::RTSolver(Solver *pS, unsigned long RTStackSize)
: RTSolverBase(pS, RTStackSize)
{
	NO_OP;
}

RTSolver::~RTSolver(void)
{
	NO_OP;
}

// write contribution to restart file
std::ostream&
RTSolver::Restart(std::ostream& out) const
{
	return out;
}

// very first setup, to be always performed
void
RTSolver::Setup(void)
{
	NO_OP;
}

// initialization to be performed only if real-time is requested
void
RTSolver::Init(void)
{
	RTSolverBase::Init();
}

// to be performed when stop is commanded by someone else
void
RTSolver::StopCommanded(void)
{
	NO_OP;
}

// write real-time related message when stop commanded by someone else
void
RTSolver::Log(void)
{
	NO_OP;
}

// wait for period to expire
void
RTSolver::Wait(void)
{
	// nanosleep
	NO_OP;
}

RTSolverBase *
ReadRTSolver(Solver *pS, MBDynParser& HP)
{
	if (HP.IsKeyWord("OS")) {
		// this
		return 0;

	}

	if (!HP.IsKeyWord("RTAI")) {
		pedantic_cerr("missing real-time type; assuming RTAI "
			"at line " << HP.GetLineData()
			<< std::endl);
	}

	return ReadRTAISolver(pS, HP);
}

/* RTSolver - end */

