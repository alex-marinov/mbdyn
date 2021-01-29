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
  * Copyright (C) 2003-2017
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * Classe che gestisce l'integrazione di un passo temporale
  *
  * Classe di base virtuale pura: NonlinearProblem 
  *
  *	metodi:
  * 		- Residual, Jacobian e Update che vengono richiesti
  * 		  da nonlinear solver
  *
  * Classe derivata: StepIntegrator
  *	 contiene:
  *		- un puntatore al data manager
  *		- le flag per il tipo di output da implementare
  *		- il numero di iterazioni
  *		- le tolleranze per la verifica della convergenza
  * 	
  *	metodi:
  *		- Predict e After Predict	
  *
  */
 
#ifndef STEPSOL_IMPL_H
#define STEPSOL_IMPL_H

#include "stepsol_tpl.h"

class CrankNicolsonIntegrator: 
	public tplStepNIntegrator<1>
{
public:
	CrankNicolsonIntegrator(const doublereal Tl, 
			const doublereal dSolTl, 
			const integer iMaxIt,
			const bool bmod_res_test);

	~CrankNicolsonIntegrator(void);

protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);

	/* Note: uses linear prediction for derivatives 
	 * (highest possible order) */
	doublereal 
	dPredDer(const doublereal dXm1mN[1],
	      const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredState(const doublereal dXm1mN[1],
		const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredDerAlg(const doublereal dXm1mN[1],
		 const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredStateAlg(const doublereal dXm1mN[1],
		   const doublereal dXP0mN[2]) const;
};


class ImplicitEulerIntegrator: 
	public tplStepNIntegrator<1>
{
public:
	ImplicitEulerIntegrator(const doublereal Tl, 
			const doublereal dSolTl, 
			const integer iMaxIt,
			const bool bmod_res_test);

	~ImplicitEulerIntegrator(void);

protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);

	/* Note: uses linear prediction for derivatives 
	 * (highest possible order) */
	doublereal 
	dPredDer(const doublereal dXm1mN[1],
	      const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredState(const doublereal dXm1mN[1],
		const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredDerAlg(const doublereal dXm1mN[1],
		 const doublereal dXP0mN[2]) const;
   
	doublereal 
	dPredStateAlg(const doublereal dXm1mN[1],
		   const doublereal dXP0mN[2]) const;
};


/* 2-step multistep (nostro metodo) - begin */

class Multistep2Solver: 
	public tplStepNIntegrator<2>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;

	doublereal m_a[2][2];
	doublereal m_b[2 + 1][2];

	doublereal m_mp[2];
	doublereal m_np[2];
 
public:
	Multistep2Solver(const doublereal Tl, 
			const doublereal dSolTol, 
			const integer iMaxIt,
			const DriveCaller* pRho,
			const DriveCaller* pAlgRho,
			const bool bmod_res_test);

	~Multistep2Solver(void);

protected:
	void SetCoef(doublereal dT, 
			doublereal dAlpha,
			enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	/* Note: uses cubic prediction for derivatives
	 * (highest possible order) */

	// dXmN: n-1, n-2, ...
	// dXP: n, n-1, n-2, ...
	doublereal 
	dPredDer(const doublereal dXm1mN[2],
			const doublereal dXP0mN[3]) const;

	// dXmN: n-1, n-2, ...
	// dXP: n, n-1, n-2, ...
	doublereal 
	dPredState(const doublereal dXm1mN[2],
			const doublereal dXP0mN[3]) const;

	// dXmN: n-1, n-2, ...
	// dXP: n, n-1, n-2, ...
	doublereal 
	dPredDerAlg(const doublereal dXm1mN[2],
			const doublereal dXP1mN[2]) const;

	// dXmN: n-1, n-2, ...
	// dXP: n, n-1, n-2, ...
	doublereal 
	dPredStateAlg(const doublereal dXm1mN[2],
			const doublereal dXP0mN[3]) const;
};

/* 2-step multistep (nostro metodo) - end */


/* Hope - begin */

class HopeSolver : 
	public tplStepNIntegrator<2>
{
protected:
	DriveOwner m_Rho;
	DriveOwner m_AlgebraicRho;
   
	bool m_bStep;
   
	doublereal m_a[2][2];
	doublereal m_b[2][2];
   
	doublereal m_mp[2];
	doublereal m_np[2];
   
public:
	HopeSolver(const doublereal Tl, 
		const doublereal dSolTol, 
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test);

	~HopeSolver(void);

protected:
	void SetCoef(doublereal dT,
		doublereal dAlpha,
		enum StepChange NewStep);

	void SetDriveHandler(const DriveHandler* pDH);

	/* Note: uses cubic prediction for derivatives
	 * (highest possible order) */
	doublereal 
	dPredDer(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;

	doublereal 
	dPredState(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;

	doublereal 
	dPredDerAlg(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;

	doublereal 
	dPredStateAlg(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const;
};

/* Hope - end */

#endif /* STEPSOL_IMPL_H */
