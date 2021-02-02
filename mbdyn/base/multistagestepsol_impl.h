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

#ifndef MULTISTAGESTEPSOL_IMPL_H
#define MULTISTAGESTEPSOL_IMPL_H

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
#include "multistagestepsol_tpl.h"

/* TunableBatheSolver2 - begin */

class TunableBatheSolver2:
	public tplStageNIntegrator<2>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;

	doublereal m_gamma;
	doublereal m_dRho;
	doublereal m_dAlgebraicRho;

	// TODO: move to tplStageNIntegrator
	doublereal m_mp[2];
	doublereal m_np[2];

public:
	TunableBatheSolver2(const doublereal Tl,
		const doublereal dSolTol,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~TunableBatheSolver2(void);

protected:
	void SetCoefForStageS(unsigned uStage,
		doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	doublereal
	dPredDerForStageS(unsigned uStage,
		const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;

	doublereal
	dPredStateForStageS(unsigned uStage,
		const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;

	doublereal
	dPredDerAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;

	doublereal
	dPredStateAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;
};

/* TunableBatheSolver2 - end */

#endif // MULTISTAGESTEPSOL_IMPL_H
