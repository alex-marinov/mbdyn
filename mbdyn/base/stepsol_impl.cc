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
  * classi che impementano l'integrazione al passo
  */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "schurdataman.h"
#include "external.h"
#include "ls.h"
#include "solver.h"
#include "invsolver.h"
#include "stepsol_impl.h"

/* CrankNicolson - begin */

CrankNicolsonIntegrator::CrankNicolsonIntegrator(const doublereal dTl,
		const doublereal dSolTl,
		const integer iMaxIt,
		const bool bmod_res_test)
: tplStepNIntegrator<1>(iMaxIt, dTl, dSolTl, bmod_res_test)
{
	NO_OP;
}

CrankNicolsonIntegrator::~CrankNicolsonIntegrator(void)
{
	NO_OP;
}

void
CrankNicolsonIntegrator::SetCoef(doublereal dT,
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
	db0Differential = db0Algebraic = dT*dAlpha/2.;
}

/* Nota: usa predizione lineare per le derivate (massimo ordine possibile) */
doublereal
CrankNicolsonIntegrator::dPredDer(const doublereal dXm1mN[1],
	const doublereal dXP0mN[2]) const
{
	return dXP0mN[IDX_XPm1];
}

doublereal
CrankNicolsonIntegrator::dPredState(const doublereal dXm1mN[1],
		const doublereal dXP0mN[2]) const
{
	return dXm1mN[IDX_Xm1] + db0Differential*(dXP0mN[IDX_XP0] + dXP0mN[IDX_XPm1]);
}

doublereal
CrankNicolsonIntegrator::dPredDerAlg(const doublereal dXm1mN[1],
		const doublereal dXP0mN[2]) const
{
	return dXP0mN[IDX_XPm1];
}

doublereal
CrankNicolsonIntegrator::dPredStateAlg(const doublereal dXm1mN[1],
		const doublereal dXP0mN[2]) const
{
	return db0Differential*(dXP0mN[IDX_XP0] + dXP0mN[IDX_XPm1]);
}

/* CrankNicolson - end */


/* Implicit Euler - begin */

ImplicitEulerIntegrator::ImplicitEulerIntegrator(const doublereal dTl,
		const doublereal dSolTl,
		const integer iMaxIt,
		const bool bmod_res_test)
: tplStepNIntegrator<1>(iMaxIt, dTl, dSolTl, bmod_res_test)
{
	NO_OP;
}

ImplicitEulerIntegrator::~ImplicitEulerIntegrator(void)
{
	NO_OP;
}

void
ImplicitEulerIntegrator::SetCoef(doublereal dT,
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
	db0Differential = db0Algebraic = dT*dAlpha;
}

/* Nota: usa predizione lineare per le derivate (massimo ordine possibile) */
doublereal
ImplicitEulerIntegrator::dPredDer(const doublereal dXm1mN[1],
	      const doublereal dXP0mN[2]) const
{
	return dXP0mN[IDX_XPm1];
}

doublereal
ImplicitEulerIntegrator::dPredState(const doublereal dXm1mN[1],
		const doublereal dXP0mN[2]) const
{
	return dXm1mN[IDX_Xm1] + db0Differential*dXP0mN[IDX_XP0];
}

doublereal
ImplicitEulerIntegrator::dPredDerAlg(const doublereal dXm1mN[1],
		const doublereal dXP0mN[2]) const
{
	return dXP0mN[IDX_XPm1];
}

doublereal
ImplicitEulerIntegrator::dPredStateAlg(const doublereal dXm1mN[1],
		const doublereal dXP0mN[2]) const
{
	return db0Differential*dXP0mN[IDX_XP0];
}

/* Implicit Euler - end */


/* 2-step multistep (nostro metodo) - begin */

Multistep2Solver::Multistep2Solver(const doublereal Tl,
		const doublereal dSolTl,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test)
: tplStepNIntegrator<2>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

Multistep2Solver::~Multistep2Solver(void)
{
	NO_OP;
}

void
Multistep2Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
Multistep2Solver::SetCoef(doublereal dT,
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
	doublereal dRho = m_Rho.dGet();
	doublereal dAlgebraicRho = m_AlgebraicRho.dGet();

	doublereal dDen = 2.*(1. + dAlpha) - (1. - dRho)*(1. - dRho);
	doublereal dBeta = dAlpha*((1. - dRho)*(1. - dRho)*(2. + dAlpha)
		+2.*(2.*dRho - 1.)*(1. + dAlpha))/dDen;
	doublereal dDelta = .5*dAlpha*dAlpha*(1. - dRho)*(1. - dRho)/dDen;

	m_mp[0] = -6.*dAlpha*dAlpha*(1. + dAlpha);
	m_mp[1] = -m_mp[0];
	m_np[0] = 1. + 4.*dAlpha + 3.*dAlpha*dAlpha;
	m_np[1] = dAlpha*(2. + 3.*dAlpha);

	m_a[IDX_A1][DIFFERENTIAL] = 1. - dBeta;
	m_a[IDX_A2][DIFFERENTIAL] = dBeta;
	m_b[IDX_B0][DIFFERENTIAL] = dT*(dDelta/dAlpha + dAlpha/2);
	m_b[IDX_B1][DIFFERENTIAL] = dT*(dBeta/2. + dAlpha/2. - dDelta/dAlpha*(1. + dAlpha));
	m_b[IDX_B2][DIFFERENTIAL] = dT*(dBeta/2. + dDelta);

	DEBUGCOUT("Predict()" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "beta  = " << dBeta << std::endl
			<< "delta = " << dDelta << std::endl
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl);

	/* Coefficienti del metodo - variabili algebriche */
	if (dAlgebraicRho != dRho) {
		dDen = 2.*(1. + dAlpha) - (1. - dAlgebraicRho)*(1. - dAlgebraicRho);
		dBeta = dAlpha*((1. - dAlgebraicRho)*(1. - dAlgebraicRho)*(2. + dAlpha)
				+2.*(2.*dAlgebraicRho - 1.)*(1. + dAlpha))/dDen;
		dDelta = .5*dAlpha*dAlpha*(1. - dAlgebraicRho)*(1. - dAlgebraicRho)/dDen;

		// m_a[IDX_A1][ALGEBRAIC] not used, since we assume the corresponding IX == 0
		m_a[IDX_A2][ALGEBRAIC] = dBeta;
		m_b[IDX_B0][ALGEBRAIC] = dT*(dDelta/dAlpha + dAlpha/2.);
		m_b[IDX_B1][ALGEBRAIC] = dT*(dBeta/2. + dAlpha/2. - dDelta/dAlpha*(1. + dAlpha));
		m_b[IDX_B2][ALGEBRAIC] = dT*(dBeta/2. + dDelta);

	} else {
		// m_a[IDX_A1][ALGEBRAIC] not used, since we assume the corresponding IX == 0
		m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
		m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
		m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
		m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
	}

	DEBUGCOUT("Algebraic coefficients:" << std::endl
			<< "beta  = " << dBeta << std::endl
			<< "delta = " << dDelta << std::endl
			<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
			<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
			<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl
			<< "b2    = " << m_b[IDX_B2][ALGEBRAIC] << std::endl);

	DEBUGCOUT("Asymptotic rho: "
			<< -m_b[IDX_B1][DIFFERENTIAL]/(2.*m_b[IDX_B0][DIFFERENTIAL]) << std::endl
			<< "Discriminant: "
			<< m_b[IDX_B1][DIFFERENTIAL]*m_b[IDX_B1][DIFFERENTIAL] - 4.*m_b[IDX_B2][DIFFERENTIAL]*m_b[IDX_B0][DIFFERENTIAL]
			<< std::endl
			<< "Asymptotic rho for algebraic variables: "
			<< -m_b[IDX_B1][ALGEBRAIC]/(2.*m_b[IDX_B0][ALGEBRAIC]) << std::endl
			<< "Discriminant: "
			<< m_b[IDX_B1][ALGEBRAIC]*m_b[IDX_B1][ALGEBRAIC] - 4.*m_b[IDX_B2][ALGEBRAIC]*m_b[IDX_B0][ALGEBRAIC]
			<< std::endl);

	/* Vengono modificati per la predizione, dopo che sono stati usati per
	 * costruire gli altri coefficienti */
	m_mp[0] /= dT;
	m_mp[1] /= dT;

	/* valori di ritorno */
	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
	//std::cout<<"PredictCoef= "<<mp[0]<<", "<<mp[1]<<", "<<np[0]<<", "<<np[1]<<std::endl;
	//std::cout<<"Coef= "<<a[0][DIFFERENTIAL]<<", "<<a[1][DIFFERENTIAL]<<", "<<b[0][DIFFERENTIAL]<<", "<<b[1][DIFFERENTIAL]<<", "<<b[2][DIFFERENTIAL]<<std::endl;
}

/* Nota: usa predizione cubica per le derivate (massimo ordine possibile) */
doublereal
Multistep2Solver::dPredDer(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const
{
	return m_mp[0]*dXm1mN[IDX_Xm1] + m_mp[1]*dXm1mN[IDX_Xm2]
		+ m_np[0]*dXP0mN[IDX_XPm1] + m_np[1]*dXP0mN[IDX_XPm2];
}

doublereal
Multistep2Solver::dPredState(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const
{
	return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1] + m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm2]
		+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0] + m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1] + m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm2];
}

// we assume that X_{n-2} = 0 and we subtract X_{n-1} (the integral of the Lagrange multipliers is inessential)
// we obtain that:
//   - coefficient n-1 is X_{n-1} - X_{n-1} = 0; 
//   - coefficient n-2 is X_{n-2} - X_{n-1} = - X_{n-1}
doublereal
Multistep2Solver::dPredDerAlg(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const
{
	return m_np[0]*dXP0mN[IDX_XPm1] + m_np[1]*dXP0mN[IDX_XPm2]
		- m_mp[1]*dXm1mN[IDX_Xm1];
}

doublereal
Multistep2Solver::dPredStateAlg(const doublereal dXm1mN[2],
		const doublereal dXP0mN[3]) const
{
	return m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0] + m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1] + m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm2]
		- m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1];
}

/* 2-step multistep (nostro metodo) - end */


/* Hope - begin */

HopeSolver::HopeSolver(const doublereal Tl,
		const doublereal dSolTl,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test)
: tplStepNIntegrator<2>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho), m_bStep(0)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

HopeSolver::~HopeSolver(void)
{
	NO_OP;
}

void
HopeSolver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
HopeSolver::SetCoef(doublereal dT,
	doublereal dAlpha,
	enum StepChange NewStep)
{
#if 0
	if (dAlpha != 1.) {
		cerr << "HOPE time step integrator is not implemented yet in variable step form" << std::endl;
		throw ErrNotImplementedYet();
	}
#endif

	if (NewStep == NEWSTEP) {
		ASSERT(m_bStep == flag(0) || m_bStep == flag(1));
		m_bStep = 1 - m_bStep;	// Commuta il valore di bStep
	}

	doublereal dTMod = dT*dAlpha;

	/* Differential coefficients */
	m_mp[0] = -6.*dAlpha*(1. + dAlpha);
	m_mp[1] = -m_mp[0];
	m_np[0] = 1. + 4.*dAlpha + 3.*dAlpha*dAlpha;
	m_np[1] = dAlpha*(2. + 3.*dAlpha);

	if (m_bStep) {
		m_b[IDX_B0][DIFFERENTIAL] = m_b[IDX_B1][DIFFERENTIAL]
			= m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B1][ALGEBRAIC]
			= db0Algebraic = db0Differential = dTMod/2.;	// dT/4.;

	} else {
		doublereal dRho = m_Rho.dGet();
		doublereal dALPHA = 4.*dRho/(3. + dRho);

		m_a[IDX_A1][DIFFERENTIAL] = (4. - dALPHA)/3.;
		m_a[IDX_A2][DIFFERENTIAL] = (dALPHA - 1.)/3.;
		m_b[IDX_B0][DIFFERENTIAL] = dTMod*(4. - dALPHA)/6.;	// dT*(4.-dALPHA)/12.;
		m_b[IDX_B1][DIFFERENTIAL] = dTMod*dALPHA/2.;		// dT*dALPHA/4.;

		DEBUGCOUT("Predict()" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "HOPE - Alpha = " << dALPHA << std::endl
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl);

		/* Coefficienti del metodo - variabili algebriche */
		doublereal dAlgebraicRho = m_AlgebraicRho.dGet();
		doublereal dAlgebraicALPHA = 4.*dAlgebraicRho/(3. + dAlgebraicRho);

		if (dAlgebraicRho != dRho) {
			// m_a[IDX_A1][ALGEBRAIC] not used, since we assume the corresponding IX == 0
			m_a[IDX_A2][ALGEBRAIC] = (dAlgebraicALPHA - 1.)/3.;
			m_b[IDX_B0][ALGEBRAIC] = dTMod*(4. - dAlgebraicALPHA)/6.; // dT*(4. - dAlgebraicALPHA)/12.;
			m_b[IDX_B1][ALGEBRAIC] = dTMod*dAlgebraicALPHA/2.; // dT*dAlgebraicALPHA/4.;

		} else {
			// m_a[IDX_A1][ALGEBRAIC] not used, since we assume the corresponding IX == 0
			m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
			m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
			m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
		}

		DEBUGCOUT("Algebraic coefficients:" << std::endl
			<< "HOPE - Alpha = " << dAlgebraicALPHA << std::endl
			<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
			<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
			<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl);

		/* valori di ritorno */
		db0Differential = m_b[IDX_B0][DIFFERENTIAL];
		db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
	}

	/* Vengono modificati per la predizione, dopo che sono stati usati per
	 * costruire gli altri coefficienti */
	m_mp[0] /= dT;
	m_mp[1] /= dT;
}

/* Nota: usa predizione cubica per le derivate (massimo ordine possibile) */
doublereal
HopeSolver::dPredDer(const doublereal dXm1mN[2],
	const doublereal dXP0mN[3]) const
{
	return m_mp[0]*dXm1mN[IDX_Xm1] + m_mp[1]*dXm1mN[IDX_Xm2]
		+ m_np[0]*dXP0mN[IDX_XPm1] + m_np[1]*dXP0mN[IDX_XPm2];
}

doublereal
HopeSolver::dPredState(const doublereal dXm1mN[2],
	const doublereal dXP0mN[3]) const
{
	if (m_bStep) {
		return dXm1mN[IDX_Xm1] + m_b[IDX_B0][ALGEBRAIC]*(dXP0mN[IDX_XP0] + dXP0mN[IDX_XPm1]);

	} else {
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1] + m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm2]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0] + m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];
	}
}

doublereal
HopeSolver::dPredDerAlg(const doublereal dXm1mN[2],
	const doublereal dXP0mN[3]) const
{
	return m_np[0]*dXP0mN[IDX_XPm1] + m_np[1]*dXP0mN[IDX_XPm2] - m_mp[1]*dXm1mN[IDX_Xm1];
}

doublereal
HopeSolver::dPredStateAlg(const doublereal dXm1mN[2],
	const doublereal dXP0mN[3]) const
{
	if (m_bStep) {
		return m_b[IDX_B0][ALGEBRAIC]*(dXP0mN[IDX_XP0] + dXP0mN[IDX_XPm1]);
	} else {
		return m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0] + m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1]
			-m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1];
	}
}

/* Hope - end */

