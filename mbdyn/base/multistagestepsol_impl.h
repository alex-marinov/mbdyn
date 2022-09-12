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

/* TunableBatheSolver - begin */

class TunableBatheSolver:
	public tplStageNIntegrator<2>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;

	doublereal m_gamma;
	doublereal m_dRho;
	doublereal m_dAlgebraicRho;

	// TODO: move to tplStageNIntegrator
	//doublereal m_mp[2];
	//doublereal m_np[2];

public:
	TunableBatheSolver(const doublereal Tl,
		const doublereal dSolTol,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~TunableBatheSolver(void);

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

/* TunableBatheSolver - end */

/* Msstc3Solver - begin */

class Msstc3Solver:
	public tplStageNIntegrator<3>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;

	doublereal m_gamma;
	doublereal m_dRho;
	doublereal m_dAlgebraicRho;

public:
	Msstc3Solver(const doublereal Tl,
		const doublereal dSolTol,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~Msstc3Solver(void);

protected:
	void SetCoefForStageS(unsigned uStage,
		doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	doublereal
	dPredDerForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;

	doublereal
	dPredStateForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;

	doublereal
	dPredDerAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;

	doublereal
	dPredStateAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;
};

/* Msstc3Solver - end */

/* Mssth3Solver - begin */

class Mssth3Solver:
	public tplStageNIntegrator<3>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;

	doublereal m_gamma;
	doublereal m_c2;
	doublereal m_dRho;
	doublereal m_dAlgebraicRho;

public:
	Mssth3Solver(const doublereal Tl,
		const doublereal dSolTol,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~Mssth3Solver(void);

protected:
	void SetCoefForStageS(unsigned uStage,
		doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	doublereal
	dPredDerForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;

	doublereal
	dPredStateForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;

	doublereal
	dPredDerAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;

	doublereal
	dPredStateAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;
};

/* Mssth3Solver - end */

/* Msstc4Solver - begin */

class Msstc4Solver:
	public tplStageNIntegrator<4>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;

	doublereal m_gamma;
	doublereal m_dRho;
	doublereal m_dAlgebraicRho;

public:
	Msstc4Solver(const doublereal Tl,
		const doublereal dSolTol,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~Msstc4Solver(void);

protected:
	void SetCoefForStageS(unsigned uStage,
		doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	doublereal
	dPredDerForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;

	doublereal
	dPredStateForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;

	doublereal
	dPredDerAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;

	doublereal
	dPredStateAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;
};

/* Msstc4Solver - end */

/* Mssth4Solver - begin */

class Mssth4Solver:
	public tplStageNIntegrator<4>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;

	doublereal m_gamma;
	doublereal m_dRho;
	doublereal m_c2;
	doublereal m_c3;
	doublereal m_dAlgebraicRho;

public:
	Mssth4Solver(const doublereal Tl,
		const doublereal dSolTol,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~Mssth4Solver(void);

protected:
	void SetCoefForStageS(unsigned uStage,
		doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	doublereal
	dPredDerForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;

	doublereal
	dPredStateForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;

	doublereal
	dPredDerAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;

	doublereal
	dPredStateAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;
};

/* Mssth4Solver - end */

/* Msstc5Solver - begin */

class Msstc5Solver:
	public tplStageNIntegrator<5>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;

	doublereal m_gamma;
	doublereal m_dRho;
	doublereal m_dAlgebraicRho;

public:
	Msstc5Solver(const doublereal Tl,
		const doublereal dSolTol,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~Msstc5Solver(void);

protected:
	void SetCoefForStageS(unsigned uStage,
		doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	doublereal
	dPredDerForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;

	doublereal
	dPredStateForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;

	doublereal
	dPredDerAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;

	doublereal
	dPredStateAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;
};

/* Msstc5Solver - end */

/* Mssth5Solver - begin */

class Mssth5Solver:
	public tplStageNIntegrator<5>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;

	doublereal m_gamma;
	doublereal m_dRho;
	doublereal m_c2;
	doublereal m_c3;
	doublereal m_c4;
	doublereal m_dAlgebraicRho;

public:
	Mssth5Solver(const doublereal Tl,
		const doublereal dSolTol,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~Mssth5Solver(void);

protected:
	void SetCoefForStageS(unsigned uStage,
		doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	doublereal
	dPredDerForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;

	doublereal
	dPredStateForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;

	doublereal
	dPredDerAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;

	doublereal
	dPredStateAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;
};

/* Mssth5Solver - end */

/*DIRK33Solver - begin */

class DIRK33Solver:
	public tplStageNIntegrator<3>
{
protected:
	//DriveOwner m_Rho;
	//DriveOwner m_AlgebraicRho;

	doublereal m_gamma;
	doublereal m_c2;
	//doublereal m_dRho;
	//doublereal m_dAlgebraicRho;

public:
	DIRK33Solver(const doublereal Tl,
		const doublereal dSolTol,
		const integer iMaxIt,
		//const DriveCaller* pRho,
		//const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~DIRK33Solver(void);

protected:
	void SetCoefForStageS(unsigned uStage,
		doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	//void SetDriveHandler(const DriveHandler* pDH);

	doublereal
	dPredDerForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;

	doublereal
	dPredStateForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;

	doublereal
	dPredDerAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;

	doublereal
	dPredStateAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const;
};

/* DIRK33Solver - end */

/* DIRK43Solver - begin */

class DIRK43Solver:
	public tplStageNIntegrator<4>
{
protected:
	//DriveOwner m_Rho;
	//DriveOwner m_AlgebraicRho;

	doublereal m_gamma;
	doublereal m_c2;
	doublereal m_c3;
	//doublereal m_dRho;
	//doublereal m_dAlgebraicRho;

public:
	DIRK43Solver(const doublereal Tl,
		const doublereal dSolTol,
		const integer iMaxIt,
		//const DriveCaller* pRho,
		//const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~DIRK43Solver(void);

protected:
	void SetCoefForStageS(unsigned uStage,
		doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	//void SetDriveHandler(const DriveHandler* pDH);

	doublereal
	dPredDerForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;

	doublereal
	dPredStateForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;

	doublereal
	dPredDerAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;

	doublereal
	dPredStateAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const;
};

/* DIRK43Solver - end */

/* DIRK54Solver - begin */

class DIRK54Solver:
	public tplStageNIntegrator<5>
{
protected:
	//DriveOwner m_Rho;
	//DriveOwner m_AlgebraicRho;

	doublereal m_gamma;
	doublereal m_c2;
	doublereal m_c3;
	doublereal m_c4;
	//doublereal m_dRho;
	//doublereal m_dAlgebraicRho;

public:
	DIRK54Solver(const doublereal Tl,
		const doublereal dSolTol,
		const integer iMaxIt,
		//const DriveCaller* pRho,
		//const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~DIRK54Solver(void);

protected:
	void SetCoefForStageS(unsigned uStage,
		doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	//void SetDriveHandler(const DriveHandler* pDH);

	doublereal
	dPredDerForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;

	doublereal
	dPredStateForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;

	doublereal
	dPredDerAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;

	doublereal
	dPredStateAlgForStageS(unsigned uStage,
		const doublereal dXm1mN[5],
		const doublereal dXP0mN[6]) const;
};

/* DIRK54Solver - end */

#endif // MULTISTAGESTEPSOL_IMPL_H
