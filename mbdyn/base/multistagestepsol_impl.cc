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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "schurdataman.h"
#include "external.h"
#include "ls.h"
#include "solver.h"
#include "invsolver.h"
#include "stepsol.h"
#include "multistagestepsol_impl.h"
#include "stepsol.hc"

/* TunableBatheSolver - begin */

TunableBatheSolver::TunableBatheSolver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	const DriveCaller* pRho,
	const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplStageNIntegrator<2>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

TunableBatheSolver::~TunableBatheSolver(void)
{
	NO_OP;
}

void
TunableBatheSolver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
TunableBatheSolver::SetCoefForStageS(unsigned uStage,
	doublereal dT,
	doublereal dAlpha,
	enum StepChange /* NewStep */)
{
	switch (uStage) {
	case 1:
		m_dRho = m_Rho.dGet();
		m_dAlgebraicRho = m_AlgebraicRho.dGet();

		if (m_dRho == 1.) {
			m_gamma = 1./2.;
		} else {
			m_gamma = (2. - sqrt(2.+ 2.*m_dRho))/(1.- m_dRho);
		}

		ASSERT(pDM != NULL);
		pDM->SetTime(pDM->dGetTime() - (1. - m_gamma)*dT, dT, pDM->pGetDrvHdl()->iGetStep());
	
		m_a[IDX_A1][DIFFERENTIAL] = 1.;
		m_a[IDX_A2][DIFFERENTIAL] = 0.;	// Unused
		m_b[IDX_B0][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B1][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B2][DIFFERENTIAL] = 0.;	// Unused

		DEBUGCOUT("PredictForStageS(1)" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = not needed" << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = not needed" << std::endl);
		break;

	case 2:
		ASSERT(pDM != NULL);
		pDM->SetTime(pDM->dGetTime() + (1. - m_gamma)*dT, dT, pDM->pGetDrvHdl()->iGetStep());
	
		/*mp[0] = -12.*dAlpha*(1.+dAlpha)/dT;
		mp[1] = -mp[0];
		np[0] = 1.+4.*dAlpha+3.*dAlpha*dAlpha;
		np[1] = dAlpha*(2.+3.*dAlpha);*/

		m_a[IDX_A1][DIFFERENTIAL] = 0.;
		m_a[IDX_A2][DIFFERENTIAL] = 1.;
		m_b[IDX_B0][DIFFERENTIAL] = (1. - m_gamma)*dT/(m_gamma*m_dRho - m_gamma + 2.);
		m_b[IDX_B1][DIFFERENTIAL] = (1. + m_dRho)*dT/(2.*(m_gamma*m_dRho - m_gamma + 2.));
		m_b[IDX_B2][DIFFERENTIAL] = (2.*m_gamma*m_dRho - m_dRho + 1.)*dT/(2.*(m_gamma*m_dRho - m_gamma + 2.));

		DEBUGCOUT("PredictForStageS(2)" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl);
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
	m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
	m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "a1    = " << m_a[IDX_A1][ALGEBRAIC] << std::endl
		<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_B2][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
}

doublereal
TunableBatheSolver::dPredDerForStageS(unsigned uStage,
	const doublereal dXm1mN[2],
	const doublereal dXP0mN[3]) const
{
	switch (uStage) {
	case 1:
		// return dXPm1;
		return dXP0mN[IDX_XPm1];

	case 2:
		// return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;
		// return dXPm1;
		return dXP0mN[IDX_XPm1]; // FIXME: at k-1 or at s1?

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
TunableBatheSolver::dPredStateForStageS(unsigned uStage,
	const doublereal dXm1mN[2],
	const doublereal dXP0mN[3]) const
{
	switch (uStage) {
	case 1:
		// return a[0][DIFFERENTIAL]*dXm1 + b[0][DIFFERENTIAL]*dXP + b[1][DIFFERENTIAL]*dXPm1;
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 2:
		// return a[0][DIFFERENTIAL]*dXm1 + a[1][DIFFERENTIAL]*dXm2
		//	+ b[0][DIFFERENTIAL]*dXP + b[1][DIFFERENTIAL]*dXPm1 + b[2][DIFFERENTIAL]*dXPm2;
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
TunableBatheSolver::dPredDerAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[2],
	const doublereal dXP0mN[3]) const
{
	switch (uStage) {
	case 1:
		// return dXPm1;
		return dXP0mN[IDX_XPm1];

	case 2:
		// return np[0]*dXPm1+np[1]*dXPm2-mp[1]*dXm1;
		// return dXPm1;
		return dXP0mN[IDX_XPm1]; // FIXME: at k-1 or at s1?

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
TunableBatheSolver::dPredStateAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[2],
	const doublereal dXP0mN[3]) const
{
	switch (uStage) {
	case 1:
		// return a[0][ALGEBRAIC]*dXm1 + b[0][ALGEBRAIC]*dXP + b[1][ALGEBRAIC]*dXPm1;
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 2:
		// return np[0]*dXPm1+np[1]*dXPm2-mp[1]*dXm1;
		// return dXPm1;
		return dXP0mN[IDX_XPm1]; // FIXME: at k-1 or s1?

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* TunableBatheSolver - end */
