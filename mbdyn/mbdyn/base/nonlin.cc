/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
 
 /* 
  *
  * Copyright (C) 2003
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */
  
#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */
 
#include <nonlin.h>  
#ifdef USE_MPI
#include <mbcomm.h>
#include <schsolman.h>
#endif /* USE_MPI */

#include <dofown.h>
#include <umfpackwrap.h>
#include <unistd.h>
#include <output.h>

SolverDiagnostics::SolverDiagnostics(unsigned OF)
{
	SetOutputFlags(OF);
}

SolverDiagnostics::~SolverDiagnostics(void)
{
	NO_OP;
}
	
void
SolverDiagnostics::SetOutputFlags(unsigned OF)
{
	OutputFlags = OF;
}

void
SolverDiagnostics::AddOutputFlags(unsigned OF)
{
	OutputFlags |= OF;
}

void
SolverDiagnostics::DelOutputFlags(unsigned OF)
{
	OutputFlags &= ~OF;
}

NonlinearSolver::NonlinearSolver(void)
: Size(0),
TotJac(0)
#ifdef USE_EXTERNAL
, ExtStepType(External::ERROR)  
#endif /* USE_EXTERNAL */
#ifdef __HACK_SCALE_RES__
, pScale(NULL) 
#endif /* __HACK_SCALE_RES__ */
{
	NO_OP;
}

#ifdef __HACK_SCALE_RES__
void
NonlinearSolver::SetScale(const VectorHandler* pScl)
{
	pScale = (VectorHandler *)pScl;
}  
#endif /* __HACK_SCALE_RES__ */

NonlinearSolver::~NonlinearSolver(void)
{
	NO_OP;
}

integer
NonlinearSolver::TotalAssembledJacobian(void)
{
	return TotJac;
}


#ifdef USE_EXTERNAL
void NonlinearSolver::SetExternal(const External::ExtMessage Ty)
{
	ExtStepType = Ty;
	return;
}

void  NonlinearSolver::SendExternal(void)
{
	switch (ExtStepType) {
		case (External::EMPTY):
			External::SendNull();
			break;

		case (External::REGULAR):
			External::SendRegular();
			break;

		case (External::CLOSE):
			External::SendClose();
			break;

		case (External::ERROR):
		default:
			External::SendError();
	}
	return;
}
#endif /* USE_EXTERNAL */
