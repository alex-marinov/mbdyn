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
 
 /* 
  *
  * Copyright (C) 2008
  * Alessandro Fumagalli	<alessandro.fumagalli@polimi.it>
  *
  * Handler of inverse dynamics problem:
  * 
  *   - Initialize data (nodes and elements)
  *   - Allocate memory for solution vectors 
  *   - Handle outputs 
  *
  */
  
#ifndef INVERSESOLVER_H
#define INVERSESOLVER_H  

#include <unistd.h>
#include <cfloat>
#include <cmath>

class InverseSolver;
#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
#include "dataman.h"
#include "schurdataman.h"
#include "schsolman.h"
#include "solver.h"
#include <deque>
#include "linsol.h"
#include "stepsol.h"
#include "nonlin.h"
#include "mfree.h"
#include "precond.h"

class InverseSolver : public Solver {
protected:
	InverseDynamics::Type ProblemType;
	doublereal dw1[3];
	doublereal dw2[3];
	MyVectorHandler* pXPrimePrime;
	MyVectorHandler* pLambda;

	bool bFullResTest;

   	/* Lettura dati */
   	void ReadData(MBDynParser& HP);

public:   
   	/* costruttore */
   	InverseSolver(MBDynParser& HP, 
		const std::string& sInputFileName, 
		const std::string& sOutputFileName,
		unsigned int nThreads,
		bool bParallel = false);

   	/* distruttore: esegue tutti i distruttori e libera la memoria */
   	virtual ~InverseSolver(void);

   	/* esegue la simulazione */
	virtual bool Prepare(void);
	virtual bool Start(void);
	virtual bool Advance(void);

   	// virtual void Run(void);

	std::ostream& Restart(std::ostream& out, DataManager::eRestart type) const;

	InverseDynamics::Type GetProblemType(void) const;
	void GetWeight(InverseDynamics::Order iOrder, doublereal& dw1, doublereal& dw2) const;
};

/* InverseSolver - end */

#endif /* INVERSESOLVER_H */

