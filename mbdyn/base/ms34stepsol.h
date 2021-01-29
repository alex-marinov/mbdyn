/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2021
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
 * Author: Huimin Zhang <huimin.zhang@polimi.it> 2021
 */

#ifndef MS34STEPSOL_H
#define MS34STEPSOL_H

#include <unistd.h>
#include <cfloat>
#include <cmath>
#include <deque>

#include "myassert.h"
#include "mynewmem.h"
#include "except.h"
#include "solman.h"
#include "dataman.h"
#include "dofown.h"
#include "drive.h"
#include "nonlinpb.h"
#include "nonlin.h"
#include "stepsol_tpl.h"

/* TunableStep3Solver - begin */

class TunableStep3Solver: 
	public tplStepNIntegrator<3>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;
   
	doublereal m_a[3][2];
	doublereal m_b[4][2];

	doublereal m_mp[2];
	doublereal m_np[2];
   
public:
	TunableStep3Solver(const doublereal Tl, 
			const doublereal dSolTol, 
			const integer iMaxIt,
			const DriveCaller* pRho,
			const DriveCaller* pAlgRho,
			const bool bmod_res_test);

	~TunableStep3Solver(void);

protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	doublereal 
	dPredDer(const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;

	doublereal 
	dPredState(const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;

	doublereal 
	dPredDerAlg(const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;
 
	doublereal 
	dPredStateAlg(const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;
};

/* TunableStep3Solver - end */


/* TunableStep4Solver - begin */

class TunableStep4Solver: 
	public tplStepNIntegrator<4>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;
   
	doublereal m_a[4][2];
	doublereal m_b[5][2];

	doublereal m_mp[2];
	doublereal m_np[2];
   
public:
	TunableStep4Solver(const doublereal Tl, 
			const doublereal dSolTol, 
			const integer iMaxIt,
			const DriveCaller* pRho,
			const DriveCaller* pAlgRho,
			const bool bmod_res_test);

	~TunableStep4Solver(void);

protected:
	void SetCoef(doublereal dT, 
		doublereal dAlpha,
		enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	doublereal 
	dPredDer(const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;
   
	doublereal 
	dPredState(const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;

	doublereal 
	dPredDerAlg(const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;

	doublereal 
	dPredStateAlg(const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;
};

/* Step4Solver - end */

#endif // MS34STEPSOL_H
