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
 
#ifndef RTSOLVER_H
#define RTSOLVER_H  

class Solver;

/* RTSolverBase - begin */

class RTSolverBase {
protected:
	Solver *pS;
	unsigned long RTStackSize;

public:
	RTSolverBase(Solver *pS, unsigned long RTStackSize);
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

/* RTSolver - begin */

class RTSolver : public RTSolverBase {
public:
	RTSolver(Solver *pS, unsigned long RTStackSize);
	~RTSolver(void);

	// write contribution to restart file
	std::ostream& Restart(std::ostream& out) const;
	// very first setup, to be always performed
	void Setup(void);
	// initialization to be performed only if real-time is requested
	void Init(void);
	// to be performed when stop is commanded by someone else
	void StopCommanded(void);
	// write real-time related message when stop commanded by someone else
	void Log(void);
	// wait for period to expire
	void Wait(void);
};

/* RTSolver - end */

extern RTSolverBase *
ReadRTSolver(Solver *pS, MBDynParser& HP);

#endif // RTSOLVER_H

