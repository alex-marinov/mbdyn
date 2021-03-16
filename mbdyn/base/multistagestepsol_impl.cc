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
			m_gamma = (2. - sqrt(2. + 2.*m_dRho))/(1. - m_dRho);
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
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (1. - m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_np[0] = 1.;
			//m_np[1] = 0.;

			//second-order prediction
			doublereal dalpha = (1. - m_gamma) / m_gamma;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((1. - m_gamma) * dT);
			m_mp[1] = -m_mp[0];
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 1.;
			m_b[IDX_B0][DIFFERENTIAL] = (1. - m_gamma) * dT / (m_gamma * m_dRho - m_gamma + 2.);
			m_b[IDX_B1][DIFFERENTIAL] = (1. + m_dRho) * dT / (2. * (m_gamma * m_dRho - m_gamma + 2.));
			m_b[IDX_B2][DIFFERENTIAL] = (2. * m_gamma * m_dRho - m_dRho + 1.) * dT / (2. * (m_gamma * m_dRho - m_gamma + 2.));

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
		}

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
		return dXP0mN[IDX_XPm1];//constant prediction 

	case 2:
		return m_mp[0]*dXm1mN[IDX_Xs1]
			+ m_mp[1]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs1]
			+ m_np[1]*dXP0mN[IDX_XPm1];

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
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 2:
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
		return dXP0mN[IDX_XPm1];

	case 2:
		return dXP0mN[IDX_XPs1]; // FIXME: at k-1 or at s1?

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
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* TunableBatheSolver - end */

/* Msstc3Solver - begin */

Msstc3Solver::Msstc3Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	const DriveCaller* pRho,
	const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplStageNIntegrator<3>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

Msstc3Solver::~Msstc3Solver(void)
{
	NO_OP;
}

void
Msstc3Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
Msstc3Solver::SetCoefForStageS(unsigned uStage,
	doublereal dT,
	doublereal dAlpha,
	enum StepChange /* NewStep */)
{
	switch (uStage) {
	case 1:
		m_dRho = m_Rho.dGet();
		m_dAlgebraicRho = m_AlgebraicRho.dGet();

		if (abs(m_dRho - 0.0) < 1.e-6)
		{
			m_gamma = 0.360850612858796;
		}else if (abs(m_dRho - 0.1) < 1.e-6)
		{
			m_gamma = 0.357238916409316;
		}else if (abs(m_dRho - 0.2) < 1.e-6)
		{
			m_gamma = 0.353891613236448;
		}else if (abs(m_dRho - 0.3) < 1.e-6)
		{
			m_gamma = 0.350771031685692;
		}else if (abs(m_dRho - 0.4) < 1.e-6)
		{
			m_gamma = 0.347847215754394;
		}else if (abs(m_dRho - 0.5) < 1.e-6)
		{
			m_gamma = 0.345095922844178;
		}else if (abs(m_dRho - 0.6) < 1.e-6)
		{
			m_gamma = 0.342497237181382;
		}else if (abs(m_dRho - 0.7) < 1.e-6)
		{
			m_gamma = 0.340034583544952;
		}else if (abs(m_dRho - 0.8) < 1.e-6)
		{
			m_gamma = 0.337694009358336;
		}else if (abs(m_dRho - 0.9) < 1.e-6)
		{
			m_gamma = 0.335463651513774;
		}else if (abs(m_dRho - 1.0) < 1.e-6)
		{
			m_gamma = 1./3.;
		}else 
		{
			silent_cerr("Please select rho in 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, and 1.0." << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ASSERT(pDM != NULL);
		pDM->SetTime(pDM->dGetTime() - (1. - m_gamma)*dT, dT, pDM->pGetDrvHdl()->iGetStep());
	
		m_a[IDX_A1][DIFFERENTIAL] = 1.;
		m_a[IDX_A2][DIFFERENTIAL] = 0.;	// Unused
		m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B0][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B1][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B2][DIFFERENTIAL] = 0.;	// Unused
		m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused

		DEBUGCOUT("PredictForStageS(1)" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = not needed" << std::endl
			<< "a3    = not needed" << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = not needed" << std::endl
			<< "b3    = not needed" << std::endl);
		break;

	case 2:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 1.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(2)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = not needed" << std::endl);
			break;
		}

	case 3:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (1. - 2. * m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;

			//second-order prediction
			doublereal dalpha = (1. - 2. * m_gamma) / m_gamma;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((1. - 2. * m_gamma) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;

			//fourth-order prediction
			//doublereal dalpha = (1. - 2.*m_gamma)/m_gamma;
			//m_mp[0] = - dalpha*dalpha*(15.*dalpha*dalpha*dalpha + 68.*dalpha*dalpha + 99.*dalpha + 46.)/(4.*(1. - 2.*m_gamma)*dT);
			//m_mp[1] = 4.*dalpha*dalpha*(dalpha*dalpha + 3.*dalpha + 2.)/((1. - 2.*m_gamma)*dT);
			//m_mp[2] = - m_mp[0] - m_mp[1];
			//m_np[0] = 5.*dalpha*dalpha*dalpha*dalpha/4. + 6.*dalpha*dalpha*dalpha + 39.*dalpha*dalpha/4. + 6.*dalpha + 1.;
			//m_np[1] = dalpha*(5.*dalpha*dalpha*dalpha + 20.*dalpha*dalpha + 24.*dalpha + 8.);
			//m_np[2] = dalpha*(5.*dalpha*dalpha*dalpha + 16.*dalpha*dalpha + 15.*dalpha + 4.)/4.;

			doublereal da1 = 1. - 3. * m_gamma / 2.;
			doublereal da2 = 1. / 2. - 3. * m_gamma / 2. + 3. * m_gamma * m_gamma / 4.;
			doublereal da3 = m_dRho * m_gamma * m_gamma * m_gamma / 8.;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 1.;
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = (m_gamma / 8. + da1 / 4. + da2 / (2. * m_gamma) + da3 / (m_gamma * m_gamma)) * dT;
			m_b[IDX_B2][DIFFERENTIAL] = (m_gamma / 2. + da1 / 2. - 2. * da3 / (m_gamma * m_gamma)) * dT;
			m_b[IDX_B3][DIFFERENTIAL] = (3. * m_gamma / 8. + da1 / 4. - da2 / (2. * m_gamma) + da3 / (m_gamma * m_gamma)) * dT;

			DEBUGCOUT("PredictForStageS(3)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl);
			break;
		}

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
	m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
	m_a[IDX_A3][ALGEBRAIC] = m_a[IDX_A3][DIFFERENTIAL];
	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
	m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
	m_b[IDX_B3][ALGEBRAIC] = m_b[IDX_B3][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "a1    = " << m_a[IDX_A1][ALGEBRAIC] << std::endl
		<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
		<< "a3    = " << m_a[IDX_A3][ALGEBRAIC] << std::endl
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_B2][ALGEBRAIC] << std::endl
		<< "b3    = " << m_b[IDX_B3][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
}

doublereal
Msstc3Solver::dPredDerForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];//constant prediction

	case 2:
		return m_mp[0]*dXm1mN[IDX_Xs1]
			+ m_mp[1]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs1]
			+ m_np[1]*dXP0mN[IDX_XPm1]; 

	case 3:
		return m_mp[0]*dXm1mN[IDX_Xs2]
			+ m_mp[1]*dXm1mN[IDX_Xs1]
			+ m_mp[2]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs2]
			+ m_np[1]*dXP0mN[IDX_XPs1]
			+ m_np[2]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Msstc3Solver::dPredStateForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Msstc3Solver::dPredDerAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return dXP0mN[IDX_XPs1]; 

	case 3:
		return dXP0mN[IDX_XPs2];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Msstc3Solver::dPredStateAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* Msstc3Solver - end */

/* Mssth3Solver - begin */

Mssth3Solver::Mssth3Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	const DriveCaller* pRho,
	const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplStageNIntegrator<3>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

Mssth3Solver::~Mssth3Solver(void)
{
	NO_OP;
}

void
Mssth3Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
Mssth3Solver::SetCoefForStageS(unsigned uStage,
	doublereal dT,
	doublereal dAlpha,
	enum StepChange /* NewStep */)
{
	switch (uStage) {
	case 1:
		m_dRho = m_Rho.dGet();
		m_dAlgebraicRho = m_AlgebraicRho.dGet();

		if (abs(m_dRho - 0.0) < 1.e-6)
		{
			m_gamma = 0.871733043016920;
		}else if (abs(m_dRho - 0.1) < 1.e-6)
		{
			m_gamma = 0.842973630818818;
		}else if (abs(m_dRho - 0.2) < 1.e-6)
		{
			m_gamma = 0.817001579025844;
		}else if (abs(m_dRho - 0.3) < 1.e-6)
		{
			m_gamma = 0.793294418242268;
		}else if (abs(m_dRho - 0.4) < 1.e-6)
		{
			m_gamma = 0.771462000921670;
		}else if (abs(m_dRho - 0.5) < 1.e-6)
		{
			m_gamma = 0.751204450030570;
		}else if (abs(m_dRho - 0.6) < 1.e-6)
		{
			m_gamma = 0.732285620206694;
		}else if (abs(m_dRho - 0.7) < 1.e-6)
		{
			m_gamma = 0.714515623934468;
		}else if (abs(m_dRho - 0.8) < 1.e-6)
		{
			m_gamma = 0.697738906149738;
		}else if (abs(m_dRho - 0.9) < 1.e-6)
		{
			m_gamma = 0.681825845543858;
		}else if (abs(m_dRho - 1.0) < 1.e-6)
		{
			m_gamma = 2./3.;
		}else 
		{
			silent_cerr("Please select rho in 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, and 1.0." << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ASSERT(pDM != NULL);
		pDM->SetTime(pDM->dGetTime() - (1. - m_gamma)*dT, dT, pDM->pGetDrvHdl()->iGetStep());
	
		m_a[IDX_A1][DIFFERENTIAL] = 1.;
		m_a[IDX_A2][DIFFERENTIAL] = 0.;	// Unused
		m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B0][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B1][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B2][DIFFERENTIAL] = 0.;	// Unused
		m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused

		DEBUGCOUT("PredictForStageS(1)" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = not needed" << std::endl
			<< "a3    = not needed" << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = not needed" << std::endl
			<< "b3    = not needed" << std::endl);
		break;

	case 2:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 1.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(2)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = not needed" << std::endl);
			break;
		}

	case 3:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (1. - 2. * m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;

			//second-order prediction
			doublereal dalpha = (1. - 2. * m_gamma) / m_gamma;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((1. - 2. * m_gamma) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;

			//fourth-order prediction
			//doublereal dalpha = (1. - 2.*m_gamma)/m_gamma;
			//m_mp[0] = - dalpha*dalpha*(15.*dalpha*dalpha*dalpha + 68.*dalpha*dalpha + 99.*dalpha + 46.)/(4.*(1. - 2.*m_gamma)*dT);
			//m_mp[1] = 4.*dalpha*dalpha*(dalpha*dalpha + 3.*dalpha + 2.)/((1. - 2.*m_gamma)*dT);
			//m_mp[2] = - m_mp[0] - m_mp[1];
			//m_np[0] = 5.*dalpha*dalpha*dalpha*dalpha/4. + 6.*dalpha*dalpha*dalpha + 39.*dalpha*dalpha/4. + 6.*dalpha + 1.;
			//m_np[1] = dalpha*(5.*dalpha*dalpha*dalpha + 20.*dalpha*dalpha + 24.*dalpha + 8.);
			//m_np[2] = dalpha*(5.*dalpha*dalpha*dalpha + 16.*dalpha*dalpha + 15.*dalpha + 4.)/4.;

			doublereal da1 = 1. - 3. * m_gamma / 2.;
			doublereal da2 = 1. / 2. - 3. * m_gamma / 2. + 3. * m_gamma * m_gamma / 4.;
			doublereal da3 = 1. / 6. - 3. * m_gamma / 4. + 3. * m_gamma * m_gamma / 4. - m_gamma * m_gamma * m_gamma / 8.;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 1.;
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = (m_gamma / 8. + da1 / 4. + da2 / (2. * m_gamma) + da3 / (m_gamma * m_gamma)) * dT;
			m_b[IDX_B2][DIFFERENTIAL] = (m_gamma / 2. + da1 / 2. - 2. * da3 / (m_gamma * m_gamma)) * dT;
			m_b[IDX_B3][DIFFERENTIAL] = (3. * m_gamma / 8. + da1 / 4. - da2 / (2. * m_gamma) + da3 / (m_gamma * m_gamma)) * dT;

			DEBUGCOUT("PredictForStageS(3)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl);
			break;
		}

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
	m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
	m_a[IDX_A3][ALGEBRAIC] = m_a[IDX_A3][DIFFERENTIAL];
	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
	m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
	m_b[IDX_B3][ALGEBRAIC] = m_b[IDX_B3][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "a1    = " << m_a[IDX_A1][ALGEBRAIC] << std::endl
		<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
		<< "a3    = " << m_a[IDX_A3][ALGEBRAIC] << std::endl
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_B2][ALGEBRAIC] << std::endl
		<< "b3    = " << m_b[IDX_B3][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
}

doublereal
Mssth3Solver::dPredDerForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return m_mp[0]*dXm1mN[IDX_Xs1]
			+ m_mp[1]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs1]
			+ m_np[1]*dXP0mN[IDX_XPm1]; 

	case 3:
		return m_mp[0]*dXm1mN[IDX_Xs2]
			+ m_mp[1]*dXm1mN[IDX_Xs1]
			+ m_mp[2]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs2]
			+ m_np[1]*dXP0mN[IDX_XPs1]
			+ m_np[2]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Mssth3Solver::dPredStateForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Mssth3Solver::dPredDerAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return dXP0mN[IDX_XPs1]; 

	case 3:
		return dXP0mN[IDX_XPs2];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Mssth3Solver::dPredStateAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm1];


	case 3:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* Mssth3Solver - end */

/* Msstc4Solver - begin */

Msstc4Solver::Msstc4Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	const DriveCaller* pRho,
	const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplStageNIntegrator<4>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

Msstc4Solver::~Msstc4Solver(void)
{
	NO_OP;
}

void
Msstc4Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
Msstc4Solver::SetCoefForStageS(unsigned uStage,
	doublereal dT,
	doublereal dAlpha,
	enum StepChange /* NewStep */)
{
	switch (uStage) {
	case 1:
		m_dRho = m_Rho.dGet();
		m_dAlgebraicRho = m_AlgebraicRho.dGet();

		if (abs(m_dRho - 0.0) < 1.e-6)
		{
			m_gamma = 0.262757473460932;
		}else if (abs(m_dRho - 0.1) < 1.e-6)
		{
			m_gamma = 0.261097241892944;
		}else if (abs(m_dRho - 0.2) < 1.e-6)
		{
			m_gamma = 0.259555167637696;
		}else if (abs(m_dRho - 0.3) < 1.e-6)
		{
			m_gamma = 0.258114414514710;
		}else if (abs(m_dRho - 0.4) < 1.e-6)
		{
			m_gamma = 0.256761609839890;
		}else if (abs(m_dRho - 0.5) < 1.e-6)
		{
			m_gamma = 0.255485941113696;
		}else if (abs(m_dRho - 0.6) < 1.e-6)
		{
			m_gamma = 0.254278531804168;
		}else if (abs(m_dRho - 0.7) < 1.e-6)
		{
			m_gamma = 0.253131998166274;
		}else if (abs(m_dRho - 0.8) < 1.e-6)
		{
			m_gamma = 0.252040127172992;
		}else if (abs(m_dRho - 0.9) < 1.e-6)
		{
			m_gamma = 0.250997637660844;
		}else if (abs(m_dRho - 1.0) < 1.e-6)
		{
			m_gamma = 1./4.;
		}else 
		{
			silent_cerr("Please select rho in 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, and 1.0." << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ASSERT(pDM != NULL);
		pDM->SetTime(pDM->dGetTime() - (1. - m_gamma)*dT, dT, pDM->pGetDrvHdl()->iGetStep());
	
		m_a[IDX_A1][DIFFERENTIAL] = 1.;
		m_a[IDX_A2][DIFFERENTIAL] = 0.;	// Unused
		m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
		m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B0][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B1][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B2][DIFFERENTIAL] = 0.;	// Unused
		m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused

		DEBUGCOUT("PredictForStageS(1)" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = not needed" << std::endl
			<< "a3    = not needed" << std::endl
			<< "a4    = not needed" << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = not needed" << std::endl
			<< "b3    = not needed" << std::endl
			<< "b4    = not needed" << std::endl);
		break;

	case 2:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;//Unused
			//m_mp[3] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;//Unused
			//m_np[3] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.; //Unused
			m_mp[3] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.; //Unused
			m_np[3] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 1.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(2)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = not needed" << std::endl
											<< "a4    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = not needed" << std::endl
											<< "b4    = not needed" << std::endl);
			break;
		}

	case 3:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 1.;
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B3][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(3)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    =  not needed" << std::endl);
			break;
		}

	case 4:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (1. - 3. * m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;

			//second-order prediction
			doublereal dalpha = (1. - 3. * m_gamma) / m_gamma;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((1. - 3. * m_gamma) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.;
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.;

			//fourth-order prediction
			//doublereal dalpha = (1.-3.*m_gamma)/m_gamma;
			//m_mp[0] = - dalpha*dalpha*(15.*dalpha*dalpha*dalpha + 68.*dalpha*dalpha + 99.*dalpha + 46.)/(4.*(1. - 3.*m_gamma)*dT);
			//m_mp[1] = 4.*dalpha*dalpha*(dalpha*dalpha + 3.*dalpha + 2.)/((1. - 3.*m_gamma)*dT);
			//m_mp[2] = - m_mp[0] - m_mp[1];
			//m_mp[3] = 0.;
			//m_np[0] = 5.*dalpha*dalpha*dalpha*dalpha/4. + 6.*dalpha*dalpha*dalpha + 39.*dalpha*dalpha/4. + 6.*dalpha + 1.;
			//m_np[1] = dalpha*(5.*dalpha*dalpha*dalpha + 20.*dalpha*dalpha + 24.*dalpha + 8.);
			//m_np[2] = dalpha*(5.*dalpha*dalpha*dalpha + 16.*dalpha*dalpha + 15.*dalpha + 4.)/4.;
			//m_np[3] = 0.;

			//sixth-order prediction
			//doublereal dalpha = (1.-3.*m_gamma)/m_gamma;
			//m_mp[0] = - dalpha*dalpha*(77.*dalpha*dalpha*dalpha*dalpha*dalpha + 774.*dalpha*dalpha*dalpha*dalpha + 3010.*dalpha*dalpha*dalpha + 5640.*dalpha*dalpha + 5073.*dalpha + 1746.)/(108.*(1. - 3.*m_gamma)*dT);
			//m_mp[1] = - dalpha*dalpha*dalpha*(7.*dalpha*dalpha*dalpha*dalpha + 60.*dalpha*dalpha*dalpha + 185.*dalpha*dalpha + 240.*dalpha + 108.)/(4.*(1. - 3.*m_gamma)*dT);
			//m_mp[2] = dalpha*dalpha*(dalpha + 3.)*(dalpha + 3.)*(7.*dalpha*dalpha*dalpha + 24.*dalpha*dalpha + 23.*dalpha + 6.)/(4.*(1. - 3.*m_gamma)*dT);
			//m_mp[3] = - m_mp[0] - m_mp[1] - m_mp[2];
			//m_np[0] = 7.*dalpha*dalpha*dalpha*dalpha*dalpha*dalpha/36. + 2.*dalpha*dalpha*dalpha*dalpha*dalpha + 145.*dalpha*dalpha*dalpha*dalpha/18. + 16.*dalpha*dalpha*dalpha + 193.*dalpha*dalpha/12. + 22.*dalpha/3. + 1.;
			//m_np[1] = dalpha*(7.*dalpha*dalpha*dalpha*dalpha*dalpha + 66.*dalpha*dalpha*dalpha*dalpha + 235.*dalpha*dalpha*dalpha + 388.*dalpha*dalpha + 288.*dalpha + 72.)/4.;
			//m_np[2] = dalpha*(7.*dalpha*dalpha*dalpha*dalpha*dalpha + 60.*dalpha*dalpha*dalpha*dalpha + 190.*dalpha*dalpha*dalpha + 272.*dalpha*dalpha + 171.*dalpha + 36.)/4.;
			//m_np[3] = dalpha*(7.*dalpha*dalpha*dalpha*dalpha*dalpha + 54.*dalpha*dalpha*dalpha*dalpha + 155.*dalpha*dalpha*dalpha + 204.*dalpha*dalpha + 120.*dalpha + 24.)/36.;

			doublereal da1 = 1. - 2. * m_gamma;
			doublereal da2 = 1. / 2. - 2. * m_gamma + 3. * m_gamma * m_gamma / 2.;
			doublereal da4 = m_dRho * m_gamma * m_gamma * m_gamma * m_gamma / 16.;
			doublereal da3 = -(3. * m_gamma * m_gamma * m_gamma * m_gamma / 8. - da2 * da2 - 2. * da4) / (2. * da1);

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.;
			m_a[IDX_A4][DIFFERENTIAL] = 1.;
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = (m_gamma / 16. + da1 / 8. + da2 / (4. * m_gamma) + da3 / (2. * m_gamma * m_gamma) + da4 / (m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B2][DIFFERENTIAL] = (5. * m_gamma / 16. + 3. * da1 / 8. + da2 / (4. * m_gamma) - da3 / (2. * m_gamma * m_gamma) - 3. * da4 / (m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B3][DIFFERENTIAL] = (11. * m_gamma / 16. + 3. * da1 / 8. - da2 / (4. * m_gamma) - da3 / (2. * m_gamma * m_gamma) + 3. * da4 / (m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B4][DIFFERENTIAL] = (7. * m_gamma / 16. + da1 / 8. - da2 / (4. * m_gamma) + da3 / (2. * m_gamma * m_gamma) - da4 / (m_gamma * m_gamma * m_gamma)) * dT;

			DEBUGCOUT("PredictForStageS(4)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = " << m_a[IDX_A4][DIFFERENTIAL] << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = " << m_b[IDX_B4][DIFFERENTIAL] << std::endl);
			break;
		}

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
	m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
	m_a[IDX_A3][ALGEBRAIC] = m_a[IDX_A3][DIFFERENTIAL];
	m_a[IDX_A4][ALGEBRAIC] = m_a[IDX_A4][DIFFERENTIAL];
	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
	m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
	m_b[IDX_B3][ALGEBRAIC] = m_b[IDX_B3][DIFFERENTIAL];
	m_b[IDX_B4][ALGEBRAIC] = m_b[IDX_B4][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "a1    = " << m_a[IDX_A1][ALGEBRAIC] << std::endl
		<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
		<< "a3    = " << m_a[IDX_A3][ALGEBRAIC] << std::endl
		<< "a4    = " << m_a[IDX_A4][ALGEBRAIC] << std::endl
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_B2][ALGEBRAIC] << std::endl
		<< "b3    = " << m_b[IDX_B3][ALGEBRAIC] << std::endl
		<< "b4    = " << m_b[IDX_B4][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
}

doublereal
Msstc4Solver::dPredDerForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];//constant prediction

	case 2:
		return m_mp[0]*dXm1mN[IDX_Xs1]
			+ m_mp[1]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs1]
			+ m_np[1]*dXP0mN[IDX_XPm1];   

	case 3:
		return m_mp[0]*dXm1mN[IDX_Xs2]
			+ m_mp[1]*dXm1mN[IDX_Xs1]
			+ m_mp[2]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs2]
			+ m_np[1]*dXP0mN[IDX_XPs1]
			+ m_np[2]*dXP0mN[IDX_XPm1]; 

	case 4:
		return m_mp[0]*dXm1mN[IDX_Xs3]
			+ m_mp[1]*dXm1mN[IDX_Xs2]
			+ m_mp[2]*dXm1mN[IDX_Xs1]
			+ m_mp[3]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs3]
			+ m_np[1]*dXP0mN[IDX_XPs2]
			+ m_np[2]*dXP0mN[IDX_XPs1]
			+ m_np[3]*dXP0mN[IDX_XPm1]; 

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Msstc4Solver::dPredStateForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Msstc4Solver::dPredDerAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return dXP0mN[IDX_XPs1]; 

	case 3:
		return dXP0mN[IDX_XPs2];

	case 4:
		return dXP0mN[IDX_XPs3];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Msstc4Solver::dPredStateAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* Msstc4Solver - end */

/* Mssth4Solver - begin */

Mssth4Solver::Mssth4Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	const DriveCaller* pRho,
	const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplStageNIntegrator<4>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

Mssth4Solver::~Mssth4Solver(void)
{
	NO_OP;
}

void
Mssth4Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
Mssth4Solver::SetCoefForStageS(unsigned uStage,
	doublereal dT,
	doublereal dAlpha,
	enum StepChange /* NewStep */)
{
	switch (uStage) {
	case 1:
		m_dRho = m_Rho.dGet();
		m_dAlgebraicRho = m_AlgebraicRho.dGet();

		if (abs(m_dRho - 0.0) < 1.e-6)
		{
			m_gamma = 1.145632124964270;
		}else if (abs(m_dRho - 0.1) < 1.e-6)
		{
			m_gamma = 1.096733289951660;
		}else if (abs(m_dRho - 0.2) < 1.e-6)
		{
			m_gamma = 1.052772913684772;
		}else if (abs(m_dRho - 0.3) < 1.e-6)
		{
			m_gamma = 1.012660237941564;
		}else if (abs(m_dRho - 0.4) < 1.e-6)
		{
			m_gamma = 0.975594949624696;
		}else if (abs(m_dRho - 0.5) < 1.e-6)
		{
			m_gamma = 0.940961155243354;
		}else if (abs(m_dRho - 0.6) < 1.e-6)
		{
			m_gamma = 0.908261570073058;
		}else if (abs(m_dRho - 0.7) < 1.e-6)
		{
			m_gamma = 0.877072379804386;
		}else if (abs(m_dRho - 0.8) < 1.e-6)
		{
			m_gamma = 0.847007532134358;
		}else if (abs(m_dRho - 0.9) < 1.e-6)
		{
			m_gamma = 0.817683732241398;
		}else if (abs(m_dRho - 1.0) < 1.e-6)
		{
			m_gamma = 0.788675134594814;
		}else 
		{
			silent_cerr("Please select rho in 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, and 1.0." << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ASSERT(pDM != NULL);
		pDM->SetTime(pDM->dGetTime() - (1. - m_gamma)*dT, dT, pDM->pGetDrvHdl()->iGetStep());
	
		m_a[IDX_A1][DIFFERENTIAL] = 1.;
		m_a[IDX_A2][DIFFERENTIAL] = 0.;	// Unused
		m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
		m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B0][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B1][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B2][DIFFERENTIAL] = 0.;	// Unused
		m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused

		DEBUGCOUT("PredictForStageS(1)" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = not needed" << std::endl
			<< "a3    = not needed" << std::endl
			<< "a4    = not needed" << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = not needed" << std::endl
			<< "b3    = not needed" << std::endl
			<< "b4    = not needed" << std::endl);
		break;

	case 2:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;//Unused
			//m_mp[3] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;//Unused
			//m_np[3] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.; //Unused
			m_mp[3] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.; //Unused
			m_np[3] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 1.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(2)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = not needed" << std::endl
											<< "a4    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = not needed" << std::endl
											<< "b4    = not needed" << std::endl);
			break;
		}

	case 3:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 1.;
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B3][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(3)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = not needed" << std::endl);
			break;
		}

	case 4:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (1. - 3. * m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;

			//second-order prediction
			//doublereal dalpha = (1. - 3.*m_gamma)/m_gamma;
			//m_mp[0] = - 6.*dalpha*dalpha*(1. + dalpha)/((1. - 3.*m_gamma)*dT);
			//m_mp[1] = - m_mp[0];
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_np[0] = 1. + 4.*dalpha + 3.*dalpha*dalpha;
			//m_np[1] = dalpha*(2. + 3.*dalpha);
			//m_np[2] = 0.;
			//m_np[3] = 0.;

			//fourth-order prediction
			doublereal dalpha = (1. - 3. * m_gamma) / m_gamma;
			m_mp[0] = -dalpha * dalpha * (15. * dalpha * dalpha * dalpha + 68. * dalpha * dalpha + 99. * dalpha + 46.) / (4. * (1. - 3. * m_gamma) * dT);
			m_mp[1] = 4. * dalpha * dalpha * (dalpha * dalpha + 3. * dalpha + 2.) / ((1. - 3. * m_gamma) * dT);
			m_mp[2] = -m_mp[0] - m_mp[1];
			m_mp[3] = 0.;
			m_np[0] = 5. * dalpha * dalpha * dalpha * dalpha / 4. + 6. * dalpha * dalpha * dalpha + 39. * dalpha * dalpha / 4. + 6. * dalpha + 1.;
			m_np[1] = dalpha * (5. * dalpha * dalpha * dalpha + 20. * dalpha * dalpha + 24. * dalpha + 8.);
			m_np[2] = dalpha * (5. * dalpha * dalpha * dalpha + 16. * dalpha * dalpha + 15. * dalpha + 4.) / 4.;
			m_np[3] = 0.;

			//sixth-order prediction
			//doublereal dalpha = (1.-3.*m_gamma)/m_gamma;
			//m_mp[0] = - dalpha*dalpha*(77.*dalpha*dalpha*dalpha*dalpha*dalpha + 774.*dalpha*dalpha*dalpha*dalpha + 3010.*dalpha*dalpha*dalpha + 5640.*dalpha*dalpha + 5073.*dalpha + 1746.)/(108.*(1. - 3.*m_gamma)*dT);
			//m_mp[1] = - dalpha*dalpha*dalpha*(7.*dalpha*dalpha*dalpha*dalpha + 60.*dalpha*dalpha*dalpha + 185.*dalpha*dalpha + 240.*dalpha + 108.)/(4.*(1. - 3.*m_gamma)*dT);
			//m_mp[2] = dalpha*dalpha*(dalpha + 3.)*(dalpha + 3.)*(7.*dalpha*dalpha*dalpha + 24.*dalpha*dalpha + 23.*dalpha + 6.)/(4.*(1. - 3.*m_gamma)*dT);
			//m_mp[3] = - m_mp[0] - m_mp[1] - m_mp[2];
			//m_np[0] = 7.*dalpha*dalpha*dalpha*dalpha*dalpha*dalpha/36. + 2.*dalpha*dalpha*dalpha*dalpha*dalpha + 145.*dalpha*dalpha*dalpha*dalpha/18. + 16.*dalpha*dalpha*dalpha + 193.*dalpha*dalpha/12. + 22.*dalpha/3. + 1.;
			//m_np[1] = dalpha*(7.*dalpha*dalpha*dalpha*dalpha*dalpha + 66.*dalpha*dalpha*dalpha*dalpha + 235.*dalpha*dalpha*dalpha + 388.*dalpha*dalpha + 288.*dalpha + 72.)/4.;
			//m_np[2] = dalpha*(7.*dalpha*dalpha*dalpha*dalpha*dalpha + 60.*dalpha*dalpha*dalpha*dalpha + 190.*dalpha*dalpha*dalpha + 272.*dalpha*dalpha + 171.*dalpha + 36.)/4.;
			//m_np[3] = dalpha*(7.*dalpha*dalpha*dalpha*dalpha*dalpha + 54.*dalpha*dalpha*dalpha*dalpha + 155.*dalpha*dalpha*dalpha + 204.*dalpha*dalpha + 120.*dalpha + 24.)/36.;

			doublereal da1 = 1. - 2. * m_gamma;
			doublereal da2 = 1. / 2. - 2. * m_gamma + 3. * m_gamma * m_gamma / 2.;
			doublereal da3 = 1. / 6. - m_gamma + 3. * m_gamma * m_gamma / 2. - m_gamma * m_gamma * m_gamma / 2.;
			doublereal da4 = 1. / 24. - m_gamma / 3. + 3. * m_gamma * m_gamma / 4. - m_gamma * m_gamma * m_gamma / 2. + m_gamma * m_gamma * m_gamma * m_gamma / 16.;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.;
			m_a[IDX_A4][DIFFERENTIAL] = 1.;
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = (m_gamma / 16. + da1 / 8. + da2 / (4. * m_gamma) + da3 / (2. * m_gamma * m_gamma) + da4 / (m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B2][DIFFERENTIAL] = (5. * m_gamma / 16. + 3. * da1 / 8. + da2 / (4. * m_gamma) - da3 / (2. * m_gamma * m_gamma) - 3. * da4 / (m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B3][DIFFERENTIAL] = (11. * m_gamma / 16. + 3. * da1 / 8. - da2 / (4. * m_gamma) - da3 / (2. * m_gamma * m_gamma) + 3. * da4 / (m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B4][DIFFERENTIAL] = (7. * m_gamma / 16. + da1 / 8. - da2 / (4. * m_gamma) + da3 / (2. * m_gamma * m_gamma) - da4 / (m_gamma * m_gamma * m_gamma)) * dT;

			DEBUGCOUT("PredictForStageS(4)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = " << m_a[IDX_A4][DIFFERENTIAL] << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = " << m_b[IDX_B4][DIFFERENTIAL] << std::endl);
			break;
		}

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
	m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
	m_a[IDX_A3][ALGEBRAIC] = m_a[IDX_A3][DIFFERENTIAL];
	m_a[IDX_A4][ALGEBRAIC] = m_a[IDX_A4][DIFFERENTIAL];
	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
	m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
	m_b[IDX_B3][ALGEBRAIC] = m_b[IDX_B3][DIFFERENTIAL];
	m_b[IDX_B4][ALGEBRAIC] = m_b[IDX_B4][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "a1    = " << m_a[IDX_A1][ALGEBRAIC] << std::endl
		<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
		<< "a3    = " << m_a[IDX_A3][ALGEBRAIC] << std::endl
		<< "a4    = " << m_a[IDX_A4][ALGEBRAIC] << std::endl
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_B2][ALGEBRAIC] << std::endl
		<< "b3    = " << m_b[IDX_B3][ALGEBRAIC] << std::endl
		<< "b4    = " << m_b[IDX_B4][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
}

doublereal
Mssth4Solver::dPredDerForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return m_mp[0]*dXm1mN[IDX_Xs1]
			+ m_mp[1]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs1]
			+ m_np[1]*dXP0mN[IDX_XPm1];

	case 3:
		return m_mp[0]*dXm1mN[IDX_Xs2]
			+ m_mp[1]*dXm1mN[IDX_Xs1]
			+ m_mp[2]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs2]
			+ m_np[1]*dXP0mN[IDX_XPs1]
			+ m_np[2]*dXP0mN[IDX_XPm1];

	case 4:
		return m_mp[0]*dXm1mN[IDX_Xs3]
			+ m_mp[1]*dXm1mN[IDX_Xs2]
			+ m_mp[2]*dXm1mN[IDX_Xs1]
			+ m_mp[3]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs3]
			+ m_np[1]*dXP0mN[IDX_XPs2]
			+ m_np[2]*dXP0mN[IDX_XPs1]
			+ m_np[3]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Mssth4Solver::dPredStateForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Mssth4Solver::dPredDerAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return dXP0mN[IDX_XPs1]; 

	case 3:
		return dXP0mN[IDX_XPs2];

	case 4:
		return dXP0mN[IDX_XPs3];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Mssth4Solver::dPredStateAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* Mssth4Solver - end */

/* Msstc5Solver - begin */

Msstc5Solver::Msstc5Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	const DriveCaller* pRho,
	const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplStageNIntegrator<5>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

Msstc5Solver::~Msstc5Solver(void)
{
	NO_OP;
}

void
Msstc5Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
Msstc5Solver::SetCoefForStageS(unsigned uStage,
	doublereal dT,
	doublereal dAlpha,
	enum StepChange /* NewStep */)
{
	switch (uStage) {
	case 1:
		m_dRho = m_Rho.dGet();
		m_dAlgebraicRho = m_AlgebraicRho.dGet();

		if (abs(m_dRho - 0.0) < 1.e-6)
		{
			m_gamma = 0.207114217840430;
		}else if (abs(m_dRho - 0.1) < 1.e-6)
		{
			m_gamma = 0.206191263023350;
		}else if (abs(m_dRho - 0.2) < 1.e-6)
		{
			m_gamma = 0.205333350050186;
		}else if (abs(m_dRho - 0.3) < 1.e-6)
		{
			m_gamma = 0.204531188984370;
		}else if (abs(m_dRho - 0.4) < 1.e-6)
		{
			m_gamma = 0.203777407759764;
		}else if (abs(m_dRho - 0.5) < 1.e-6)
		{
			m_gamma = 0.203066050295748;
		}else if (abs(m_dRho - 0.6) < 1.e-6)
		{
			m_gamma = 0.202392230146362;
		}else if (abs(m_dRho - 0.7) < 1.e-6)
		{
			m_gamma = 0.201751884891614;
		}else if (abs(m_dRho - 0.8) < 1.e-6)
		{
			m_gamma = 0.201141597837490;
		}else if (abs(m_dRho - 0.9) < 1.e-6)
		{
			m_gamma = 0.200558465909484;
		}else if (abs(m_dRho - 1.0) < 1.e-6)
		{
			m_gamma = 1./5.;
		}else 
		{
			silent_cerr("Please select rho in 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, and 1.0." << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ASSERT(pDM != NULL);
		pDM->SetTime(pDM->dGetTime() - (1. - m_gamma)*dT, dT, pDM->pGetDrvHdl()->iGetStep());
	
		m_a[IDX_A1][DIFFERENTIAL] = 1.;
		m_a[IDX_A2][DIFFERENTIAL] = 0.;	// Unused
		m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
		m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
		m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B0][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B1][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B2][DIFFERENTIAL] = 0.;	// Unused
		m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

		DEBUGCOUT("PredictForStageS(1)" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = not needed" << std::endl
			<< "a3    = not needed" << std::endl
			<< "a4    = not needed" << std::endl
			<< "a5    = not needed" << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = not needed" << std::endl
			<< "b3    = not needed" << std::endl
			<< "b4    = not needed" << std::endl
			<< "b5    = not needed" << std::endl);
		break;

	case 2:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;//Unused
			//m_mp[3] = 0.;//Unused
			//m_mp[4] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;//Unused
			//m_np[3] = 0.;//Unused
			//m_np[4] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.; //Unused
			m_mp[3] = 0.; //Unused
			m_mp[4] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.; //Unused
			m_np[3] = 0.; //Unused
			m_np[4] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 1.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(2)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = not needed" << std::endl
											<< "a4    = not needed" << std::endl
											<< "a5    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = not needed" << std::endl
											<< "b4    = not needed" << std::endl
											<< "b5    = not needed" << std::endl);
			break;
		}

	case 3:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;//Unused
			//m_mp[4] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;//Unused
			//m_np[4] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.; //Unused
			m_mp[4] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.; //Unused
			m_np[4] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 1.;
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B3][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(3)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = not needed" << std::endl
											<< "a5    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = not needed" << std::endl
											<< "b5    = not needed" << std::endl);
			break;
		}

	case 4:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_np[4] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;
			//m_np[4] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.;
			m_mp[4] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.;
			m_np[4] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.;
			m_a[IDX_A4][DIFFERENTIAL] = 1.;
			m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B3][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B4][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(4)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = " << m_a[IDX_A4][DIFFERENTIAL] << std::endl
											<< "a5    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = " << m_b[IDX_B4][DIFFERENTIAL] << std::endl
											<< "b5    = not needed" << std::endl);
			break;
		}

	case 5:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (1. - 4. * m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_np[4] = 0.;
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;
			//m_np[4] = 0.;

			//second-order prediction
			doublereal dalpha = (1. - 4. * m_gamma) / m_gamma;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((1. - 4. * m_gamma) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.;
			m_mp[4] = 0.;
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.;
			m_np[4] = 0.;

			//fourth-order prediction
			//doublereal dalpha = (1. - 4.*m_gamma)/m_gamma;
			//m_mp[0] = - dalpha*dalpha*(15.*dalpha*dalpha*dalpha + 68.*dalpha*dalpha + 99.*dalpha + 46.)/(4.*(1. - 4.*m_gamma)*dT);
			//m_mp[1] = 4.*dalpha*dalpha*(dalpha*dalpha + 3.*dalpha + 2.)/((1. - 4.*m_gamma)*dT);
			//m_mp[2] = - m_mp[0] - m_mp[1];
			//m_mp[3] = 0.;
			//m_mp[4] = 0.;
			//m_np[0] = 5.*dalpha*dalpha*dalpha*dalpha/4. + 6.*dalpha*dalpha*dalpha + 39.*dalpha*dalpha/4. + 6.*dalpha + 1.;
			//m_np[1] = dalpha*(5.*dalpha*dalpha*dalpha + 20.*dalpha*dalpha + 24.*dalpha + 8.);
			//m_np[2] = dalpha*(5.*dalpha*dalpha*dalpha + 16.*dalpha*dalpha + 15.*dalpha + 4.)/4.;
			//m_np[3] = 0.;
			//m_np[4] = 0.;

			//sixth-order prediction
			//doublereal dalpha = (1. - 4.*m_gamma)/m_gamma;
			//m_mp[0] = - dalpha*dalpha*(77.*dalpha*dalpha*dalpha*dalpha*dalpha + 774.*dalpha*dalpha*dalpha*dalpha + 3010.*dalpha*dalpha*dalpha + 5640.*dalpha*dalpha + 5073.*dalpha + 1746.)/(108.*(1. - 4.*m_gamma)*dT);
			//m_mp[1] = - dalpha*dalpha*dalpha*(7.*dalpha*dalpha*dalpha*dalpha + 60.*dalpha*dalpha*dalpha + 185.*dalpha*dalpha + 240.*dalpha + 108.)/(4.*(1. - 4.*m_gamma)*dT);
			//m_mp[2] = dalpha*dalpha*(dalpha + 3.)*(dalpha + 3.)*(7.*dalpha*dalpha*dalpha + 24.*dalpha*dalpha + 23.*dalpha + 6.)/(4.*(1. - 4.*m_gamma)*dT);
			//m_mp[3] = - m_mp[0] - m_mp[1] - m_mp[2];
			//m_mp[4] = 0.;
			//m_np[0] = 7.*dalpha*dalpha*dalpha*dalpha*dalpha*dalpha/36. + 2.*dalpha*dalpha*dalpha*dalpha*dalpha + 145.*dalpha*dalpha*dalpha*dalpha/18. + 16.*dalpha*dalpha*dalpha + 193.*dalpha*dalpha/12. + 22.*dalpha/3. + 1.;
			//m_np[1] = dalpha*(7.*dalpha*dalpha*dalpha*dalpha*dalpha + 66.*dalpha*dalpha*dalpha*dalpha + 235.*dalpha*dalpha*dalpha + 388.*dalpha*dalpha + 288.*dalpha + 72.)/4.;
			//m_np[2] = dalpha*(7.*dalpha*dalpha*dalpha*dalpha*dalpha + 60.*dalpha*dalpha*dalpha*dalpha + 190.*dalpha*dalpha*dalpha + 272.*dalpha*dalpha + 171.*dalpha + 36.)/4.;
			//m_np[3] = dalpha*(7.*dalpha*dalpha*dalpha*dalpha*dalpha + 54.*dalpha*dalpha*dalpha*dalpha + 155.*dalpha*dalpha*dalpha + 204.*dalpha*dalpha + 120.*dalpha + 24.)/36.;
			//m_np[4] = 0.;

			doublereal da1 = 1. - 5. * m_gamma / 2.;
			doublereal da2 = 5. * m_gamma * m_gamma / 2. - 5. * m_gamma / 2. + 1. / 2.;
			doublereal da3 = -(50. * m_gamma * m_gamma * m_gamma - 70. * m_gamma * m_gamma + 30. * m_gamma - 4. +
							   sqrt(2. * ((5. * m_dRho + 805.) * m_gamma * m_gamma * m_gamma * m_gamma * m_gamma * m_gamma +
										  (-2. * m_dRho - 2050.) * m_gamma * m_gamma * m_gamma * m_gamma * m_gamma + 2160. * m_gamma * m_gamma * m_gamma * m_gamma -
										  1200. * m_gamma * m_gamma * m_gamma + 370. * m_gamma * m_gamma - 60. * m_gamma + 4.))) /
							 8.;
			doublereal da4 = -(45. * m_gamma * m_gamma * m_gamma * m_gamma - 100. * m_gamma * m_gamma * m_gamma + 70. * m_gamma * m_gamma + (40. * da3 - 20.) * m_gamma - 16. * da3 + 2.) / 16.;
			doublereal da5 = (m_dRho * m_gamma * m_gamma * m_gamma * m_gamma * m_gamma) / 32.;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.;
			m_a[IDX_A4][DIFFERENTIAL] = 0.;
			m_a[IDX_A5][DIFFERENTIAL] = 1.;
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = (m_gamma / 32. + da1 / 16. + da2 / (8. * m_gamma) + da3 / (4. * m_gamma * m_gamma) + da4 / (2. * m_gamma * m_gamma * m_gamma) + da5 / (m_gamma * m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B2][DIFFERENTIAL] = (3. * m_gamma / 16. + da1 / 4. + da2 / (4. * m_gamma) - da4 / (m_gamma * m_gamma * m_gamma) - 4. * da5 / (m_gamma * m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B3][DIFFERENTIAL] = (m_gamma / 2. + 3. * da1 / 8. - da3 / (2. * m_gamma * m_gamma) + 6. * da5 / (m_gamma * m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B4][DIFFERENTIAL] = (13. * m_gamma / 16. + da1 / 4. - da2 / (4. * m_gamma) + da4 / (m_gamma * m_gamma * m_gamma) - 4. * da5 / (m_gamma * m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B5][DIFFERENTIAL] = (15. * m_gamma / 32. + da1 / 16. - da2 / (8. * m_gamma) + da3 / (4. * m_gamma * m_gamma) - da4 / (2. * m_gamma * m_gamma * m_gamma) + da5 / (m_gamma * m_gamma * m_gamma * m_gamma)) * dT;

			DEBUGCOUT("PredictForStageS(5)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = " << m_a[IDX_A4][DIFFERENTIAL] << std::endl
											<< "a5    = " << m_a[IDX_A5][DIFFERENTIAL] << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = " << m_b[IDX_B4][DIFFERENTIAL] << std::endl
											<< "b5    = " << m_b[IDX_B5][DIFFERENTIAL] << std::endl);
			break;
		}

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
	m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
	m_a[IDX_A3][ALGEBRAIC] = m_a[IDX_A3][DIFFERENTIAL];
	m_a[IDX_A4][ALGEBRAIC] = m_a[IDX_A4][DIFFERENTIAL];
	m_a[IDX_A5][ALGEBRAIC] = m_a[IDX_A5][DIFFERENTIAL];
	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
	m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
	m_b[IDX_B3][ALGEBRAIC] = m_b[IDX_B3][DIFFERENTIAL];
	m_b[IDX_B4][ALGEBRAIC] = m_b[IDX_B4][DIFFERENTIAL];
	m_b[IDX_B5][ALGEBRAIC] = m_b[IDX_B5][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "a1    = " << m_a[IDX_A1][ALGEBRAIC] << std::endl
		<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
		<< "a3    = " << m_a[IDX_A3][ALGEBRAIC] << std::endl
		<< "a4    = " << m_a[IDX_A4][ALGEBRAIC] << std::endl
		<< "a5    = " << m_a[IDX_A5][ALGEBRAIC] << std::endl
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_B2][ALGEBRAIC] << std::endl
		<< "b3    = " << m_b[IDX_B3][ALGEBRAIC] << std::endl
		<< "b4    = " << m_b[IDX_B4][ALGEBRAIC] << std::endl
		<< "b5    = " << m_b[IDX_B5][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
}

doublereal
Msstc5Solver::dPredDerForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return m_mp[0]*dXm1mN[IDX_Xs1]
			+ m_mp[1]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs1]
			+ m_np[1]*dXP0mN[IDX_XPm1]; 

	case 3:
		return m_mp[0]*dXm1mN[IDX_Xs2]
			+ m_mp[1]*dXm1mN[IDX_Xs1]
			+ m_mp[2]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs2]
			+ m_np[1]*dXP0mN[IDX_XPs1]
			+ m_np[2]*dXP0mN[IDX_XPm1];

	case 4:
		return m_mp[0]*dXm1mN[IDX_Xs3]
			+ m_mp[1]*dXm1mN[IDX_Xs2]
			+ m_mp[2]*dXm1mN[IDX_Xs1]
			+ m_mp[3]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs3]
			+ m_np[1]*dXP0mN[IDX_XPs2]
			+ m_np[2]*dXP0mN[IDX_XPs1]
			+ m_np[3]*dXP0mN[IDX_XPm1];

	case 5:
		return m_mp[0]*dXm1mN[IDX_Xs4]
			+ m_mp[1]*dXm1mN[IDX_Xs3]
			+ m_mp[2]*dXm1mN[IDX_Xs2]
			+ m_mp[3]*dXm1mN[IDX_Xs1]
			+ m_mp[4]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs4]
			+ m_np[1]*dXP0mN[IDX_XPs3]
			+ m_np[2]*dXP0mN[IDX_XPs2]
			+ m_np[3]*dXP0mN[IDX_XPs1]
			+ m_np[4]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Msstc5Solver::dPredStateForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 5:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs4]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A4][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A5][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B4][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B5][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Msstc5Solver::dPredDerAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return dXP0mN[IDX_XPs1]; 

	case 3:
		return dXP0mN[IDX_XPs2];

	case 4:
		return dXP0mN[IDX_XPs3];

	case 5:
		return dXP0mN[IDX_XPs4];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Msstc5Solver::dPredStateAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 5:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs4]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A4][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A5][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B4][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B5][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* Msstc5Solver - end */

/* Mssth5Solver - begin */

Mssth5Solver::Mssth5Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	const DriveCaller* pRho,
	const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplStageNIntegrator<5>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

Mssth5Solver::~Mssth5Solver(void)
{
	NO_OP;
}

void
Mssth5Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
Mssth5Solver::SetCoefForStageS(unsigned uStage,
	doublereal dT,
	doublereal dAlpha,
	enum StepChange /* NewStep */)
{
	switch (uStage) {
	case 1:
		m_dRho = m_Rho.dGet();
		m_dAlgebraicRho = m_AlgebraicRho.dGet();

		if (abs(m_dRho - 0.0) < 1.e-6)
		{
			m_gamma = 0.556107682272900;
		}else if (abs(m_dRho - 0.1) < 1.e-6)
		{
			m_gamma = 0.548282612063736;
		}else if (abs(m_dRho - 0.2) < 1.e-6)
		{
			m_gamma = 0.540919773549164;
		}else if (abs(m_dRho - 0.3) < 1.e-6)
		{
			m_gamma = 0.533956087851302;
		}else if (abs(m_dRho - 0.4) < 1.e-6)
		{
			m_gamma = 0.527340463423212;
		}else if (abs(m_dRho - 0.5) < 1.e-6)
		{
			m_gamma = 0.521030833214110;
		}else if (abs(m_dRho - 0.6) < 1.e-6)
		{
			m_gamma = 0.514992059713350;
		}else if (abs(m_dRho - 0.7) < 1.e-6)
		{
			m_gamma = 0.509194416340266;
		}else if (abs(m_dRho - 0.8) < 1.e-6)
		{
			m_gamma = 0.503612462367700;
		}else if (abs(m_dRho - 0.9) < 1.e-6)
		{
			m_gamma = 0.498224193059260;
		}else if (abs(m_dRho - 1.0) < 1.e-6)
		{
			m_gamma = 0.493010386285640;
		}else 
		{
			silent_cerr("Please select rho in 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, and 1.0." << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		

		ASSERT(pDM != NULL);
		pDM->SetTime(pDM->dGetTime() - (1. - m_gamma)*dT, dT, pDM->pGetDrvHdl()->iGetStep());
	
		m_a[IDX_A1][DIFFERENTIAL] = 1.;
		m_a[IDX_A2][DIFFERENTIAL] = 0.;	// Unused
		m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
		m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
		m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B0][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B1][DIFFERENTIAL] = m_gamma*dT/2.;
		m_b[IDX_B2][DIFFERENTIAL] = 0.;	// Unused
		m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

		DEBUGCOUT("PredictForStageS(1)" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = not needed" << std::endl
			<< "a3    = not needed" << std::endl
			<< "a4    = not needed" << std::endl
			<< "a5    = not needed" << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = not needed" << std::endl
			<< "b3    = not needed" << std::endl
			<< "b4    = not needed" << std::endl
			<< "b5    = not needed" << std::endl);
		break;

	case 2:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;//Unused
			//m_mp[3] = 0.;//Unused
			//m_mp[4] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;//Unused
			//m_np[3] = 0.;//Unused
			//m_np[4] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.; //Unused
			m_mp[3] = 0.; //Unused
			m_mp[4] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.; //Unused
			m_np[3] = 0.; //Unused
			m_np[4] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 1.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(2)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = not needed" << std::endl
											<< "a4    = not needed" << std::endl
											<< "a5    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = not needed" << std::endl
											<< "b4    = not needed" << std::endl
											<< "b5    = not needed" << std::endl);
			break;
		}

	case 3:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;//Unused
			//m_mp[4] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;//Unused
			//m_np[4] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.; //Unused
			m_mp[4] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.; //Unused
			m_np[4] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 1.;
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B3][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(3)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = not needed" << std::endl
											<< "a5    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = not needed" << std::endl
											<< "b5    = not needed" << std::endl);
			break;
		}

	case 4:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + m_gamma * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_np[4] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;
			//m_np[4] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = 1.;
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / (m_gamma * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.;
			m_mp[4] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.;
			m_np[4] = 0.; //Unused

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.;
			m_a[IDX_A4][DIFFERENTIAL] = 1.;
			m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B2][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B3][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B4][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(4)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = " << m_a[IDX_A4][DIFFERENTIAL] << std::endl
											<< "a5    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = " << m_b[IDX_B4][DIFFERENTIAL] << std::endl
											<< "b5    = not needed" << std::endl);
			break;
		}

	case 5:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (1. - 4. * m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_np[4] = 0.;
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;
			//m_np[4] = 0.;

			//second-order prediction
			//doublereal dalpha = (1. - 4.*m_gamma)/m_gamma;
			//m_mp[0] = -6.*dalpha*dalpha*(1. + dalpha)/((1. - 4.*m_gamma)*dT);
			//m_mp[1] = -m_mp[0];
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_mp[4] = 0.;
			//m_np[0] = 1. + 4.*dalpha + 3.*dalpha*dalpha;
			//m_np[1] = dalpha*(2. + 3.*dalpha);
			//m_np[2] = 0.;
			//m_np[3] = 0.;
			//m_np[4] = 0.;

			//fourth-order prediction
			//doublereal dalpha = (1. - 4.*m_gamma)/m_gamma;
			//m_mp[0] = - dalpha*dalpha*(15.*dalpha*dalpha*dalpha + 68.*dalpha*dalpha + 99.*dalpha + 46.)/(4.*(1. - 4.*m_gamma)*dT);
			//m_mp[1] = 4.*dalpha*dalpha*(dalpha*dalpha + 3.*dalpha + 2.)/((1. - 4.*m_gamma)*dT);
			//m_mp[2] = - m_mp[0] - m_mp[1];
			//m_mp[3] = 0.;
			//m_mp[4] = 0.;
			//m_np[0] = 5.*dalpha*dalpha*dalpha*dalpha/4. + 6.*dalpha*dalpha*dalpha + 39.*dalpha*dalpha/4. + 6.*dalpha + 1.;
			//m_np[1] = dalpha*(5.*dalpha*dalpha*dalpha + 20.*dalpha*dalpha + 24.*dalpha + 8.);
			//m_np[2] = dalpha*(5.*dalpha*dalpha*dalpha + 16.*dalpha*dalpha + 15.*dalpha + 4.)/4.;
			//m_np[3] = 0.;
			//m_np[4] = 0.;

			//sixth-order prediction
			doublereal dalpha = (1. - 4. * m_gamma) / m_gamma;
			m_mp[0] = -dalpha * dalpha * (77. * dalpha * dalpha * dalpha * dalpha * dalpha + 774. * dalpha * dalpha * dalpha * dalpha + 3010. * dalpha * dalpha * dalpha + 5640. * dalpha * dalpha + 5073. * dalpha + 1746.) / (108. * (1. - 4. * m_gamma) * dT);
			m_mp[1] = -dalpha * dalpha * dalpha * (7. * dalpha * dalpha * dalpha * dalpha + 60. * dalpha * dalpha * dalpha + 185. * dalpha * dalpha + 240. * dalpha + 108.) / (4. * (1. - 4. * m_gamma) * dT);
			m_mp[2] = dalpha * dalpha * (dalpha + 3.) * (dalpha + 3.) * (7. * dalpha * dalpha * dalpha + 24. * dalpha * dalpha + 23. * dalpha + 6.) / (4. * (1. - 4. * m_gamma) * dT);
			m_mp[3] = -m_mp[0] - m_mp[1] - m_mp[2];
			m_mp[4] = 0.;
			m_np[0] = 7. * dalpha * dalpha * dalpha * dalpha * dalpha * dalpha / 36. + 2. * dalpha * dalpha * dalpha * dalpha * dalpha + 145. * dalpha * dalpha * dalpha * dalpha / 18. + 16. * dalpha * dalpha * dalpha + 193. * dalpha * dalpha / 12. + 22. * dalpha / 3. + 1.;
			m_np[1] = dalpha * (7. * dalpha * dalpha * dalpha * dalpha * dalpha + 66. * dalpha * dalpha * dalpha * dalpha + 235. * dalpha * dalpha * dalpha + 388. * dalpha * dalpha + 288. * dalpha + 72.) / 4.;
			m_np[2] = dalpha * (7. * dalpha * dalpha * dalpha * dalpha * dalpha + 60. * dalpha * dalpha * dalpha * dalpha + 190. * dalpha * dalpha * dalpha + 272. * dalpha * dalpha + 171. * dalpha + 36.) / 4.;
			m_np[3] = dalpha * (7. * dalpha * dalpha * dalpha * dalpha * dalpha + 54. * dalpha * dalpha * dalpha * dalpha + 155. * dalpha * dalpha * dalpha + 204. * dalpha * dalpha + 120. * dalpha + 24.) / 36.;
			m_np[4] = 0.;

			doublereal da1 = 1. - 5. * m_gamma / 2.;
			doublereal da2 = 1. / 2. - 5. * m_gamma / 2. + 5. * m_gamma * m_gamma / 2.;
			doublereal da3 = 1. / 6. - 5. * m_gamma / 4. + 5. * m_gamma * m_gamma / 2. - 5. * m_gamma * m_gamma * m_gamma / 4.;
			doublereal da4 = 1. / 24. - 5. * m_gamma / 12. + 5. * m_gamma * m_gamma / 4. - 5. * m_gamma * m_gamma * m_gamma / 4. + 5. * m_gamma * m_gamma * m_gamma * m_gamma / 16.;
			doublereal da5 = 1. / 120. - 5. * m_gamma / 48. + 5. * m_gamma * m_gamma / 12. - 5. * m_gamma * m_gamma * m_gamma / 8. + 5. * m_gamma * m_gamma * m_gamma * m_gamma / 16. - m_gamma * m_gamma * m_gamma * m_gamma * m_gamma / 32.;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.;
			m_a[IDX_A4][DIFFERENTIAL] = 0.;
			m_a[IDX_A5][DIFFERENTIAL] = 1.;
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT / 2.;
			m_b[IDX_B1][DIFFERENTIAL] = (m_gamma / 32. + da1 / 16. + da2 / (8. * m_gamma) + da3 / (4. * m_gamma * m_gamma) + da4 / (2. * m_gamma * m_gamma * m_gamma) + da5 / (m_gamma * m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B2][DIFFERENTIAL] = (3. * m_gamma / 16. + da1 / 4. + da2 / (4. * m_gamma) - da4 / (m_gamma * m_gamma * m_gamma) - 4. * da5 / (m_gamma * m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B3][DIFFERENTIAL] = (m_gamma / 2. + 3. * da1 / 8. - da3 / (2. * m_gamma * m_gamma) + 6. * da5 / (m_gamma * m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B4][DIFFERENTIAL] = (13. * m_gamma / 16. + da1 / 4. - da2 / (4. * m_gamma) + da4 / (m_gamma * m_gamma * m_gamma) - 4. * da5 / (m_gamma * m_gamma * m_gamma * m_gamma)) * dT;
			m_b[IDX_B5][DIFFERENTIAL] = (15. * m_gamma / 32. + da1 / 16. - da2 / (8. * m_gamma) + da3 / (4. * m_gamma * m_gamma) - da4 / (2. * m_gamma * m_gamma * m_gamma) + da5 / (m_gamma * m_gamma * m_gamma * m_gamma)) * dT;

			DEBUGCOUT("PredictForStageS(5)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = " << m_a[IDX_A4][DIFFERENTIAL] << std::endl
											<< "a5    = " << m_a[IDX_A5][DIFFERENTIAL] << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = " << m_b[IDX_B4][DIFFERENTIAL] << std::endl
											<< "b5    = " << m_b[IDX_B5][DIFFERENTIAL] << std::endl);
			break;
		}

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
	m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
	m_a[IDX_A3][ALGEBRAIC] = m_a[IDX_A3][DIFFERENTIAL];
	m_a[IDX_A4][ALGEBRAIC] = m_a[IDX_A4][DIFFERENTIAL];
	m_a[IDX_A5][ALGEBRAIC] = m_a[IDX_A5][DIFFERENTIAL];
	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
	m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
	m_b[IDX_B3][ALGEBRAIC] = m_b[IDX_B3][DIFFERENTIAL];
	m_b[IDX_B4][ALGEBRAIC] = m_b[IDX_B4][DIFFERENTIAL];
	m_b[IDX_B5][ALGEBRAIC] = m_b[IDX_B5][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "a1    = " << m_a[IDX_A1][ALGEBRAIC] << std::endl
		<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
		<< "a3    = " << m_a[IDX_A3][ALGEBRAIC] << std::endl
		<< "a4    = " << m_a[IDX_A4][ALGEBRAIC] << std::endl
		<< "a5    = " << m_a[IDX_A5][ALGEBRAIC] << std::endl
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_B2][ALGEBRAIC] << std::endl
		<< "b3    = " << m_b[IDX_B3][ALGEBRAIC] << std::endl
		<< "b4    = " << m_b[IDX_B4][ALGEBRAIC] << std::endl
		<< "b5    = " << m_b[IDX_B5][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
}

doublereal
Mssth5Solver::dPredDerForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return m_mp[0]*dXm1mN[IDX_Xs1]
			+ m_mp[1]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs1]
			+ m_np[1]*dXP0mN[IDX_XPm1]; 

	case 3:
		return m_mp[0]*dXm1mN[IDX_Xs2]
			+ m_mp[1]*dXm1mN[IDX_Xs1]
			+ m_mp[2]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs2]
			+ m_np[1]*dXP0mN[IDX_XPs1]
			+ m_np[2]*dXP0mN[IDX_XPm1];

	case 4:
		return m_mp[0]*dXm1mN[IDX_Xs3]
			+ m_mp[1]*dXm1mN[IDX_Xs2]
			+ m_mp[2]*dXm1mN[IDX_Xs1]
			+ m_mp[3]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs3]
			+ m_np[1]*dXP0mN[IDX_XPs2]
			+ m_np[2]*dXP0mN[IDX_XPs1]
			+ m_np[3]*dXP0mN[IDX_XPm1];

	case 5:
		return m_mp[0]*dXm1mN[IDX_Xs4]
			+ m_mp[1]*dXm1mN[IDX_Xs3]
			+ m_mp[2]*dXm1mN[IDX_Xs2]
			+ m_mp[3]*dXm1mN[IDX_Xs1]
			+ m_mp[4]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs4]
			+ m_np[1]*dXP0mN[IDX_XPs3]
			+ m_np[2]*dXP0mN[IDX_XPs2]
			+ m_np[3]*dXP0mN[IDX_XPs1]
			+ m_np[4]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Mssth5Solver::dPredStateForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 5:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs4]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A4][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A5][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B4][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B5][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Mssth5Solver::dPredDerAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return dXP0mN[IDX_XPs1]; 

	case 3:
		return dXP0mN[IDX_XPs2];

	case 4:
		return dXP0mN[IDX_XPs3];

	case 5:
		return dXP0mN[IDX_XPs4];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
Mssth5Solver::dPredStateAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 5:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs4]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A4][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A5][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B4][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B5][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* Mssth5Solver - end */

/* DIRK33Solver - begin */

DIRK33Solver::DIRK33Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	//const DriveCaller* pRho,
	//const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplStageNIntegrator<3>(iMaxIt, Tl, dSolTl, bmod_res_test)
{
	NO_OP;
}

DIRK33Solver::~DIRK33Solver(void)
{
	NO_OP;
}


void
DIRK33Solver::SetCoefForStageS(unsigned uStage,
	doublereal dT,
	doublereal dAlpha,
	enum StepChange /* NewStep */)
{
	switch (uStage) {
	case 1:
		m_gamma = 0.43586652150845899941601945;

		ASSERT(pDM != NULL);
		pDM->SetTime(pDM->dGetTime() - (1. - 2.*m_gamma)*dT, dT, pDM->pGetDrvHdl()->iGetStep());
	
		m_a[IDX_A1][DIFFERENTIAL] = 1.;
		m_a[IDX_A2][DIFFERENTIAL] = 0.;	// Unused
		m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
		m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
		m_b[IDX_B2][DIFFERENTIAL] = 0.;	// Unused
		m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused

		DEBUGCOUT("PredictForStageS(1)" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = not needed" << std::endl
			<< "a3    = not needed" << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = not needed" << std::endl
			<< "b3    = not needed" << std::endl);
		break;

	case 2:
		{
			m_c2 = 3. / 5.;
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (m_c2 - 2. * m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = (m_c2 - 2. * m_gamma) / (2. * m_gamma);
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((m_c2 - 2. * m_gamma) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.; //Unused

			doublereal a21 = m_c2 * (m_c2 - 2. * m_gamma) / (4. * m_gamma);
			doublereal a20 = m_c2 - m_gamma - a21;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 1.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B1][DIFFERENTIAL] = a21 * dT;
			m_b[IDX_B2][DIFFERENTIAL] = a20 * dT;
			m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(2)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = not needed" << std::endl);
			break;
		}

	case 3:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (1. - m_c2) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;

			//second-order prediction
			doublereal dalpha = (1. - m_c2) / (m_c2 - 2. * m_gamma);
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((1. - m_c2) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;

			//fourth-order prediction
			//doublereal dalpha1 = (1. - m_c2) / (m_c2 - 2. * m_gamma);
			//doublereal dalpha2 = (1. - m_c2) / (2. * m_gamma);
			//m_mp[0] = -(2. * dalpha1 * dalpha1 * (1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
			//		+ dalpha1 * dalpha1 * (10. * dalpha2 * dalpha2 * dalpha2 + 25. * dalpha2 * dalpha2 + 13. * dalpha2) 
			//		+ dalpha1 * (20. * dalpha2 * dalpha2 * dalpha2 + 20. * dalpha2 * dalpha2) + 10. * dalpha2 * dalpha2 * dalpha2)) 
			//		/ ((dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (1. - m_c2) * dT);
			//m_mp[1] = (2. * (1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
			//		+ dalpha1 * dalpha1 * (-5. * dalpha2 * dalpha2 * dalpha2 + dalpha2 * dalpha2 + 4. * dalpha2) 
			//		+ dalpha1 * (-7. * dalpha2 * dalpha2 * dalpha2 - dalpha2 * dalpha2) - 2. * dalpha2 * dalpha2 * dalpha2)) 
			//		/ (dalpha1 * (1. - m_c2) * dT);
			//m_mp[2] = (2. * dalpha2 * dalpha2 * dalpha2 * dalpha2 * (1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (10. * dalpha2 + 10.) 
			//		+ dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 20. * dalpha2 + 5.) + dalpha1 * (7. * dalpha2 * dalpha2 + 7. * dalpha2) 
			//		+ 2. * dalpha2 * dalpha2)) / (dalpha1 * (dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (1. - m_c2) * dT);
			//m_np[0] = ((1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
			//		+ dalpha1 * dalpha1 * (11. * dalpha2 * dalpha2 + 10. * dalpha2 + 1.) + dalpha1 * (7. * dalpha2 * dalpha2 + 2. * dalpha2) 
			//		+ dalpha2 * dalpha2)) / ((dalpha1 + dalpha2) * (dalpha1 + dalpha2));
			//m_np[1] = (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
			//		+ dalpha1 * dalpha1 * (12. * dalpha2 * dalpha2 + 12. * dalpha2 + 2.) + dalpha1 * (9. * dalpha2 * dalpha2 + 4. * dalpha2) 
			//		+ 2. * dalpha2 * dalpha2) / dalpha1;
			//m_np[2] = (dalpha2 * dalpha2 * dalpha2 * (1. + dalpha1) * (dalpha1 * dalpha1 * (5. * dalpha2 + 4.) 
			//		+ dalpha1 * (7. * dalpha2 + 2.) + 2. * dalpha2)) / (dalpha1 * (dalpha1 + dalpha2) * (dalpha1 + dalpha2));

			doublereal q1 = (-2. + 3. * m_c2 + 6. * m_gamma * (1. - m_c2)) / (12. * m_gamma * (m_c2 - 2. * m_gamma));
			doublereal q2 = (1. - 6. * m_gamma + 6. * m_gamma * m_gamma) / (3. * m_c2 * (m_c2 - 2. * m_gamma));
			doublereal q0 = 1. - m_gamma - q1 - q2;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 1.;
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B1][DIFFERENTIAL] = q2 * dT;
			m_b[IDX_B2][DIFFERENTIAL] = q1 * dT;
			m_b[IDX_B3][DIFFERENTIAL] = q0 * dT;

			DEBUGCOUT("PredictForStageS(3)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl);
			break;
		}

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
	m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
	m_a[IDX_A3][ALGEBRAIC] = m_a[IDX_A3][DIFFERENTIAL];
	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
	m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
	m_b[IDX_B3][ALGEBRAIC] = m_b[IDX_B3][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "a1    = " << m_a[IDX_A1][ALGEBRAIC] << std::endl
		<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
		<< "a3    = " << m_a[IDX_A3][ALGEBRAIC] << std::endl
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_B2][ALGEBRAIC] << std::endl
		<< "b3    = " << m_b[IDX_B3][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
}

doublereal
DIRK33Solver::dPredDerForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];//constant prediction

	case 2:
		return m_mp[0]*dXm1mN[IDX_Xs1]
			+ m_mp[1]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs1]
			+ m_np[1]*dXP0mN[IDX_XPm1]; 

	case 3:
		return m_mp[0]*dXm1mN[IDX_Xs2]
			+ m_mp[1]*dXm1mN[IDX_Xs1]
			+ m_mp[2]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs2]
			+ m_np[1]*dXP0mN[IDX_XPs1]
			+ m_np[2]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
DIRK33Solver::dPredStateForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
DIRK33Solver::dPredDerAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return dXP0mN[IDX_XPs1]; 

	case 3:
		return dXP0mN[IDX_XPs2];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
DIRK33Solver::dPredStateAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[3],
	const doublereal dXP0mN[4]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* DIRK33Solver - end */

/* DIRK43Solver - begin */

DIRK43Solver::DIRK43Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	//const DriveCaller* pRho,
	//const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplStageNIntegrator<4>(iMaxIt, Tl, dSolTl, bmod_res_test)
{
	NO_OP;
}

DIRK43Solver::~DIRK43Solver(void)
{
	NO_OP;
}

void
DIRK43Solver::SetCoefForStageS(unsigned uStage,
	doublereal dT,
	doublereal dAlpha,
	enum StepChange /* NewStep */)
{
	switch (uStage) {
	case 1:
		m_gamma = 9. / 40.;

		ASSERT(pDM != NULL);
		pDM->SetTime(pDM->dGetTime() - (1. - 2. * m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

		m_a[IDX_A1][DIFFERENTIAL] = 1.;
		m_a[IDX_A2][DIFFERENTIAL] = 0.;	// Unused
		m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
		m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
		m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
		m_b[IDX_B2][DIFFERENTIAL] = 0.;	// Unused
		m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused

		DEBUGCOUT("PredictForStageS(1)" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = not needed" << std::endl
			<< "a3    = not needed" << std::endl
			<< "a4    = not needed" << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = not needed" << std::endl
			<< "b3    = not needed" << std::endl
			<< "b4    = not needed" << std::endl);
		break;

	case 2:
		{
			m_c2 = 9. * (2. + sqrt(2.)) / 40.;
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (m_c2 - 2. * m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;//Unused
			//m_mp[3] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;//Unused
			//m_np[3] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = (m_c2 - 2. * m_gamma) / (2. * m_gamma);
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((m_c2 - 2. * m_gamma) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.; //Unused
			m_mp[3] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.; //Unused
			m_np[3] = 0.; //Unused

			doublereal a21 = (m_c2 * (m_c2 - 2. * m_gamma)) / (4. * m_gamma);
			doublereal a20 = m_c2 - m_gamma - a21;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 1.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B1][DIFFERENTIAL] = a21 * dT;
			m_b[IDX_B2][DIFFERENTIAL] = a20 * dT;
			m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(2)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = not needed" << std::endl
											<< "a4    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = not needed" << std::endl
											<< "b4    = not needed" << std::endl);
			break;
		}

	case 3:
		{
			m_c3 = 3. / 5.;
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (m_c3 - m_c2) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = (m_c3 - m_c2) / (m_c2 - 2. * m_gamma);
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((m_c3 - m_c2) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.; //Unused

			doublereal a31 = (-4. * m_gamma * m_gamma * m_gamma + 12. * m_c3 * m_gamma * m_gamma - 4. * m_c3 * m_c3 * m_gamma 
							- 2. * m_c2 * m_c3 * m_gamma + m_c2 * m_c3 * m_c3) / (4. * m_gamma * (m_c2 - 2. * m_gamma));
			doublereal a32 = m_gamma * (2. * m_gamma * m_gamma - 4. * m_c3 * m_gamma + m_c3 * m_c3) / (m_c2 * (m_c2 - 2. * m_gamma));
			doublereal a30 = m_c3 - m_gamma - a31 - a32;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 1.;
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B1][DIFFERENTIAL] = a32 * dT;
			m_b[IDX_B2][DIFFERENTIAL] = a31 * dT;
			m_b[IDX_B3][DIFFERENTIAL] = a30 * dT;
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(3)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    =  not needed" << std::endl);
			break;
		}

	case 4:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (1. - m_c3) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;

			//second-order prediction
			doublereal dalpha = (1. - m_c3) / (m_c3 - m_c2);
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((1. - m_c3) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.;
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.;

			//fourth-order prediction
			//doublereal dalpha1 = (1. - m_c3) / (m_c3 - m_c2);
			//doublereal dalpha2 = (1. - m_c3) / (m_c2 - 2. * m_gamma);
			//m_mp[0] = -(2. * dalpha1 * dalpha1 * (1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
			//		+ dalpha1 * dalpha1 * (10. * dalpha2 * dalpha2 * dalpha2 + 25. * dalpha2 * dalpha2 + 13. * dalpha2) 
			//		+ dalpha1 * (20. * dalpha2 * dalpha2 * dalpha2 + 20. * dalpha2 * dalpha2) + 10. * dalpha2 * dalpha2 * dalpha2)) 
			//		/ ((dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (1. - m_c3) * dT);
			//m_mp[1] = (2. * (1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
			//		+ dalpha1 * dalpha1 * (-5. * dalpha2 * dalpha2 * dalpha2 + dalpha2 * dalpha2 + 4. * dalpha2) 
			//		+ dalpha1 * (-7. * dalpha2 * dalpha2 * dalpha2 - dalpha2 * dalpha2) - 2. * dalpha2 * dalpha2 * dalpha2)) 
			//		/ (dalpha1 * (1. - m_c3) * dT);
			//m_mp[2] = (2. * dalpha2 * dalpha2 * dalpha2 * dalpha2 * (1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (10. * dalpha2 + 10.) 
			//		+ dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 20. * dalpha2 + 5.) + dalpha1 * (7. * dalpha2 * dalpha2 + 7. * dalpha2) 
			//		+ 2. * dalpha2 * dalpha2)) / (dalpha1 * (dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (1. - m_c3) * dT);
			//m_mp[3] = 0.;
			//m_np[0] = ((1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
			//		+ dalpha1 * dalpha1 * (11. * dalpha2 * dalpha2 + 10. * dalpha2 + 1.) + dalpha1 * (7. * dalpha2 * dalpha2 + 2. * dalpha2) 
			//		+ dalpha2 * dalpha2)) / ((dalpha1 + dalpha2) * (dalpha1 + dalpha2));
			//m_np[1] = (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
			//		+ dalpha1 * dalpha1 * (12. * dalpha2 * dalpha2 + 12. * dalpha2 + 2.) + dalpha1 * (9. * dalpha2 * dalpha2 + 4. * dalpha2) 
			//		+ 2. * dalpha2 * dalpha2) / dalpha1;
			//m_np[2] = (dalpha2 * dalpha2 * dalpha2 * (1. + dalpha1) * (dalpha1 * dalpha1 * (5. * dalpha2 + 4.) 
			//		+ dalpha1 * (7. * dalpha2 + 2.) + 2. * dalpha2)) / (dalpha1 * (dalpha1 + dalpha2) * (dalpha1 + dalpha2));
			//m_np[3] = 0.;

			doublereal q1 = -(12. * (m_c3 - 1.) * (m_c3 - m_c2 + 1.) * m_gamma * m_gamma * m_gamma + 2. * (12. * m_c3 - 3. * m_c2 + 6. * m_c2 * m_c3 - 18. * m_c3 * m_c3 + 2.) 
							* m_gamma * m_gamma + 2. * m_c3 * (6. * m_c3 - 3. * m_c2 + 3. * m_c2 * m_c3 - 4.) * m_gamma - m_c2 * m_c3 * (3. * m_c3 - 2.)) 
							/ (12. * m_gamma * (m_c2 - 2. * m_gamma) * (m_c3 * m_c3 - 4. * m_c3 * m_gamma + 2. * m_gamma * m_gamma));
			doublereal q2 = (m_gamma * (12. * (1. - m_c3) * m_gamma * m_gamma * m_gamma + 6. * (m_c3 * m_c3 + 2. * m_c3 - 2.) * m_gamma * m_gamma - 
							2. * (6. * m_c3 * m_c3 - 3. * m_c3 - 1.) * m_gamma + m_c3 * (3. * m_c3 - 2.))) 
							/ (3. * m_c2 * (m_c2 - 2. * m_gamma) * (m_c3 * m_c3 - 4. * m_c3 * m_gamma + 2. * m_gamma * m_gamma));
			doublereal q3 = -(6. * m_gamma * m_gamma * m_gamma - 18. * m_gamma * m_gamma + 9. * m_gamma - 1.) / (3. * (m_c3 * m_c3 - 4. * m_c3 * m_gamma + 2. * m_gamma * m_gamma));
			doublereal q0 = 1. - m_gamma - q1 - q2 - q3;
			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.;
			m_a[IDX_A4][DIFFERENTIAL] = 1.;
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B1][DIFFERENTIAL] = q3 * dT;
			m_b[IDX_B2][DIFFERENTIAL] = q2 * dT;
			m_b[IDX_B3][DIFFERENTIAL] = q1 * dT;
			m_b[IDX_B4][DIFFERENTIAL] = q0 * dT;

			DEBUGCOUT("PredictForStageS(4)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = " << m_a[IDX_A4][DIFFERENTIAL] << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = " << m_b[IDX_B4][DIFFERENTIAL] << std::endl);
			break;
		}

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
	m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
	m_a[IDX_A3][ALGEBRAIC] = m_a[IDX_A3][DIFFERENTIAL];
	m_a[IDX_A4][ALGEBRAIC] = m_a[IDX_A4][DIFFERENTIAL];
	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
	m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
	m_b[IDX_B3][ALGEBRAIC] = m_b[IDX_B3][DIFFERENTIAL];
	m_b[IDX_B4][ALGEBRAIC] = m_b[IDX_B4][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "a1    = " << m_a[IDX_A1][ALGEBRAIC] << std::endl
		<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
		<< "a3    = " << m_a[IDX_A3][ALGEBRAIC] << std::endl
		<< "a4    = " << m_a[IDX_A4][ALGEBRAIC] << std::endl
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_B2][ALGEBRAIC] << std::endl
		<< "b3    = " << m_b[IDX_B3][ALGEBRAIC] << std::endl
		<< "b4    = " << m_b[IDX_B4][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
}

doublereal
DIRK43Solver::dPredDerForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];//constant prediction

	case 2:
		return m_mp[0]*dXm1mN[IDX_Xs1]
			+ m_mp[1]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs1]
			+ m_np[1]*dXP0mN[IDX_XPm1];   

	case 3:
		return m_mp[0]*dXm1mN[IDX_Xs2]
			+ m_mp[1]*dXm1mN[IDX_Xs1]
			+ m_mp[2]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs2]
			+ m_np[1]*dXP0mN[IDX_XPs1]
			+ m_np[2]*dXP0mN[IDX_XPm1]; 

	case 4:
		return m_mp[0]*dXm1mN[IDX_Xs3]
			+ m_mp[1]*dXm1mN[IDX_Xs2]
			+ m_mp[2]*dXm1mN[IDX_Xs1]
			+ m_mp[3]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs3]
			+ m_np[1]*dXP0mN[IDX_XPs2]
			+ m_np[2]*dXP0mN[IDX_XPs1]
			+ m_np[3]*dXP0mN[IDX_XPm1]; 

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
DIRK43Solver::dPredStateForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
DIRK43Solver::dPredDerAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return dXP0mN[IDX_XPs1]; 

	case 3:
		return dXP0mN[IDX_XPs2];

	case 4:
		return dXP0mN[IDX_XPs3];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
DIRK43Solver::dPredStateAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[4],
	const doublereal dXP0mN[5]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* DIRK43Solver - end */

/* DIRK54Solver - begin */

DIRK54Solver::DIRK54Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	//const DriveCaller* pRho,
	//const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplStageNIntegrator<5>(iMaxIt, Tl, dSolTl, bmod_res_test)
{
	NO_OP;
}

DIRK54Solver::~DIRK54Solver(void)
{
	NO_OP;
}


void
DIRK54Solver::SetCoefForStageS(unsigned uStage,
	doublereal dT,
	doublereal dAlpha,
	enum StepChange /* NewStep */)
{
	switch (uStage) {
	case 1:
		m_gamma = 31. / 125.;

		ASSERT(pDM != NULL);
		pDM->SetTime(pDM->dGetTime() - (1. - 2. * m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

		m_a[IDX_A1][DIFFERENTIAL] = 1.;
		m_a[IDX_A2][DIFFERENTIAL] = 0.;	// Unused
		m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
		m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
		m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
		m_b[IDX_B1][DIFFERENTIAL] = m_gamma * dT;
		m_b[IDX_B2][DIFFERENTIAL] = 0.;	// Unused
		m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused
		m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

		DEBUGCOUT("PredictForStageS(1)" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
			<< "a2    = not needed" << std::endl
			<< "a3    = not needed" << std::endl
			<< "a4    = not needed" << std::endl
			<< "a5    = not needed" << std::endl
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
			<< "b2    = not needed" << std::endl
			<< "b3    = not needed" << std::endl
			<< "b4    = not needed" << std::endl
			<< "b5    = not needed" << std::endl);
		break;

	case 2:
		{
			m_c2 = 486119545908. / 3346201505189.;
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (m_c2 - 2. * m_gamma) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;//Unused
			//m_mp[3] = 0.;//Unused
			//m_mp[4] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;//Unused
			//m_np[3] = 0.;//Unused
			//m_np[4] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = (m_c2 - 2. * m_gamma) / (2. * m_gamma);
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((m_c2 - 2. * m_gamma) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.; //Unused
			m_mp[3] = 0.; //Unused
			m_mp[4] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.; //Unused
			m_np[3] = 0.; //Unused
			m_np[4] = 0.; //Unused

			doublereal a21 = -360286518617. / 7014585480527.;
			doublereal a20 = m_c2 - m_gamma - a21;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 1.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B1][DIFFERENTIAL] = a21 * dT;
			m_b[IDX_B2][DIFFERENTIAL] = a20 * dT;
			m_b[IDX_B3][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(2)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = not needed" << std::endl
											<< "a4    = not needed" << std::endl
											<< "a5    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = not needed" << std::endl
											<< "b4    = not needed" << std::endl
											<< "b5    = not needed" << std::endl);
			break;
		}

	case 3:
		{
			m_c3 = 1043. / 1706.;
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (m_c3 - m_c2) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;//Unused
			//m_mp[4] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;//Unused
			//m_np[4] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = (m_c3 - m_c2) / (m_c2 - 2. * m_gamma);
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((m_c3 - m_c2) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.; //Unused
			m_mp[4] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.; //Unused
			m_np[4] = 0.; //Unused

			doublereal a31 = -506388693497. / 5937754990171.;
			doublereal a32 = 7149918333491. / 13390931526268.;
			doublereal a30 = m_c3 - m_gamma - a31 - a32;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 1.;
			m_a[IDX_A4][DIFFERENTIAL] = 0.; // Unused
			m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B1][DIFFERENTIAL] = a32 * dT;
			m_b[IDX_B2][DIFFERENTIAL] = a31 * dT;
			m_b[IDX_B3][DIFFERENTIAL] = a30 * dT;
			m_b[IDX_B4][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(3)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = not needed" << std::endl
											<< "a5    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = not needed" << std::endl
											<< "b5    = not needed" << std::endl);
			break;
		}

	case 4:
		{
			m_c4 = 1361. / 1300.;
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (m_c4 - m_c3) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_np[4] = 0.;//Unused
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;
			//m_np[4] = 0.;//Unused

			//second-order prediction
			doublereal dalpha = (m_c4 - m_c3) / (m_c3 - m_c2);
			m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((m_c4 - m_c3) * dT);
			m_mp[1] = -m_mp[0];
			m_mp[2] = 0.;
			m_mp[3] = 0.;
			m_mp[4] = 0.; //Unused
			m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			m_np[1] = dalpha * (2. + 3. * dalpha);
			m_np[2] = 0.;
			m_np[3] = 0.;
			m_np[4] = 0.; //Unused

			doublereal a41 = -7628305438933. / 11061539393788.;
			doublereal a42 = 21592626537567. / 14352247503901.;
			doublereal a43 = 11630056083252. / 17263101053231.;
			doublereal a40 = m_c4 - m_gamma - a41 - a42 - a43;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.;
			m_a[IDX_A4][DIFFERENTIAL] = 1.;
			m_a[IDX_A5][DIFFERENTIAL] = 0.; // Unused
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B1][DIFFERENTIAL] = a43 * dT;
			m_b[IDX_B2][DIFFERENTIAL] = a42 * dT;
			m_b[IDX_B3][DIFFERENTIAL] = a41 * dT;
			m_b[IDX_B4][DIFFERENTIAL] = a40 * dT;
			m_b[IDX_B5][DIFFERENTIAL] = 0.; // Unused

			DEBUGCOUT("PredictForStageS(4)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = " << m_a[IDX_A4][DIFFERENTIAL] << std::endl
											<< "a5    = not needed" << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = " << m_b[IDX_B4][DIFFERENTIAL] << std::endl
											<< "b5    = not needed" << std::endl);
			break;
		}

	case 5:
		{
			ASSERT(pDM != NULL);
			pDM->SetTime(pDM->dGetTime() + (1. - m_c4) * dT, dT, pDM->pGetDrvHdl()->iGetStep());

			//constant prediction
			//m_mp[0] = 0.;
			//m_mp[1] = 0.;
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_np[4] = 0.;
			//m_np[0] = 1.;
			//m_np[1] = 0.;
			//m_np[2] = 0.;
			//m_np[3] = 0.;
			//m_np[4] = 0.;

			//second-order prediction
			//doublereal dalpha = (1. - m_c4) / (m_c4 - m_c3);
			//m_mp[0] = -6. * dalpha * dalpha * (1. + dalpha) / ((1. - m_c4) * dT);
			//m_mp[1] = -m_mp[0];
			//m_mp[2] = 0.;
			//m_mp[3] = 0.;
			//m_mp[4] = 0.;
			//m_np[0] = 1. + 4. * dalpha + 3. * dalpha * dalpha;
			//m_np[1] = dalpha * (2. + 3. * dalpha);
			//m_np[2] = 0.;
			//m_np[3] = 0.;
			//m_np[4] = 0.;

			//fourth-order prediction
			doublereal dalpha1 = (1. - m_c4) / (m_c4 - m_c3);
			doublereal dalpha2 = (1. - m_c4) / (m_c3 - m_c2);
			m_mp[0] = -(2. * dalpha1 * dalpha1 * (1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
					+ dalpha1 * dalpha1 * (10. * dalpha2 * dalpha2 * dalpha2 + 25. * dalpha2 * dalpha2 + 13. * dalpha2) 
					+ dalpha1 * (20. * dalpha2 * dalpha2 * dalpha2 + 20. * dalpha2 * dalpha2) + 10. * dalpha2 * dalpha2 * dalpha2)) 
					/ ((dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (1. - m_c4) * dT);
			m_mp[1] = (2. * (1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
					+ dalpha1 * dalpha1 * (-5. * dalpha2 * dalpha2 * dalpha2 + dalpha2 * dalpha2 + 4. * dalpha2) 
					+ dalpha1 * (-7. * dalpha2 * dalpha2 * dalpha2 - dalpha2 * dalpha2) - 2. * dalpha2 * dalpha2 * dalpha2)) 
					/ (dalpha1 * (1. - m_c4) * dT);
			m_mp[2] = (2. * dalpha2 * dalpha2 * dalpha2 * dalpha2 * (1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (10. * dalpha2 + 10.) 
					+ dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 20. * dalpha2 + 5.) + dalpha1 * (7. * dalpha2 * dalpha2 + 7. * dalpha2) 
					+ 2. * dalpha2 * dalpha2)) / (dalpha1 * (dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (dalpha1 + dalpha2) * (1. - m_c4) * dT);
			m_mp[3] = 0.;
			m_mp[4] = 0.;
			m_np[0] = ((1. + dalpha1) * (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
					+ dalpha1 * dalpha1 * (11. * dalpha2 * dalpha2 + 10. * dalpha2 + 1.) + dalpha1 * (7. * dalpha2 * dalpha2 + 2. * dalpha2) 
					+ dalpha2 * dalpha2)) / ((dalpha1 + dalpha2) * (dalpha1 + dalpha2));
			m_np[1] = (dalpha1 * dalpha1 * dalpha1 * (5. * dalpha2 * dalpha2 + 8. * dalpha2 + 3.) 
					+ dalpha1 * dalpha1 * (12. * dalpha2 * dalpha2 + 12. * dalpha2 + 2.) + dalpha1 * (9. * dalpha2 * dalpha2 + 4. * dalpha2) 
					+ 2. * dalpha2 * dalpha2) / dalpha1;
			m_np[2] = (dalpha2 * dalpha2 * dalpha2 * (1. + dalpha1) * (dalpha1 * dalpha1 * (5. * dalpha2 + 4.) 
					+ dalpha1 * (7. * dalpha2 + 2.) + 2. * dalpha2)) / (dalpha1 * (dalpha1 + dalpha2) * (dalpha1 + dalpha2));
			m_np[3] = 0.;
			m_np[4] = 0.;

			doublereal q1 = -12917657251. / 5222094901039.;
			doublereal q2 = 5602338284630. / 15643096342197.;
			doublereal q3 = 9002339615474. / 18125249312447.;
			doublereal q4 = -2420307481369. / 24731958684496.;
			doublereal q0 = 1. - m_gamma - q1 - q2 - q3 - q4;

			m_a[IDX_A1][DIFFERENTIAL] = 0.;
			m_a[IDX_A2][DIFFERENTIAL] = 0.;
			m_a[IDX_A3][DIFFERENTIAL] = 0.;
			m_a[IDX_A4][DIFFERENTIAL] = 0.;
			m_a[IDX_A5][DIFFERENTIAL] = 1.;
			m_b[IDX_B0][DIFFERENTIAL] = m_gamma * dT;
			m_b[IDX_B1][DIFFERENTIAL] = q4 * dT;
			m_b[IDX_B2][DIFFERENTIAL] = q3 * dT;
			m_b[IDX_B3][DIFFERENTIAL] = q2 * dT;
			m_b[IDX_B4][DIFFERENTIAL] = q1 * dT;
			m_b[IDX_B5][DIFFERENTIAL] = q0 * dT;

			DEBUGCOUT("PredictForStageS(5)" << std::endl
											<< "Alpha = " << dAlpha << std::endl
											<< "Differential coefficients:" << std::endl
											<< "Asymptotic rho =" << m_dRho << std::endl
											<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
											<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
											<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
											<< "a4    = " << m_a[IDX_A4][DIFFERENTIAL] << std::endl
											<< "a5    = " << m_a[IDX_A5][DIFFERENTIAL] << std::endl
											<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
											<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
											<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
											<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
											<< "b4    = " << m_b[IDX_B4][DIFFERENTIAL] << std::endl
											<< "b5    = " << m_b[IDX_B5][DIFFERENTIAL] << std::endl);
			break;
		}

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
	m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
	m_a[IDX_A3][ALGEBRAIC] = m_a[IDX_A3][DIFFERENTIAL];
	m_a[IDX_A4][ALGEBRAIC] = m_a[IDX_A4][DIFFERENTIAL];
	m_a[IDX_A5][ALGEBRAIC] = m_a[IDX_A5][DIFFERENTIAL];
	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
	m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
	m_b[IDX_B3][ALGEBRAIC] = m_b[IDX_B3][DIFFERENTIAL];
	m_b[IDX_B4][ALGEBRAIC] = m_b[IDX_B4][DIFFERENTIAL];
	m_b[IDX_B5][ALGEBRAIC] = m_b[IDX_B5][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "a1    = " << m_a[IDX_A1][ALGEBRAIC] << std::endl
		<< "a2    = " << m_a[IDX_A2][ALGEBRAIC] << std::endl
		<< "a3    = " << m_a[IDX_A3][ALGEBRAIC] << std::endl
		<< "a4    = " << m_a[IDX_A4][ALGEBRAIC] << std::endl
		<< "a5    = " << m_a[IDX_A5][ALGEBRAIC] << std::endl
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_B1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_B2][ALGEBRAIC] << std::endl
		<< "b3    = " << m_b[IDX_B3][ALGEBRAIC] << std::endl
		<< "b4    = " << m_b[IDX_B4][ALGEBRAIC] << std::endl
		<< "b5    = " << m_b[IDX_B5][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];
}

doublereal
DIRK54Solver::dPredDerForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return m_mp[0]*dXm1mN[IDX_Xs1]
			+ m_mp[1]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs1]
			+ m_np[1]*dXP0mN[IDX_XPm1]; 

	case 3:
		return m_mp[0]*dXm1mN[IDX_Xs2]
			+ m_mp[1]*dXm1mN[IDX_Xs1]
			+ m_mp[2]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs2]
			+ m_np[1]*dXP0mN[IDX_XPs1]
			+ m_np[2]*dXP0mN[IDX_XPm1];

	case 4:
		return m_mp[0]*dXm1mN[IDX_Xs3]
			+ m_mp[1]*dXm1mN[IDX_Xs2]
			+ m_mp[2]*dXm1mN[IDX_Xs1]
			+ m_mp[3]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs3]
			+ m_np[1]*dXP0mN[IDX_XPs2]
			+ m_np[2]*dXP0mN[IDX_XPs1]
			+ m_np[3]*dXP0mN[IDX_XPm1];

	case 5:
		return m_mp[0]*dXm1mN[IDX_Xs4]
			+ m_mp[1]*dXm1mN[IDX_Xs3]
			+ m_mp[2]*dXm1mN[IDX_Xs2]
			+ m_mp[3]*dXm1mN[IDX_Xs1]
			+ m_mp[4]*dXm1mN[IDX_Xm1]
			+ m_np[0]*dXP0mN[IDX_XPs4]
			+ m_np[1]*dXP0mN[IDX_XPs3]
			+ m_np[2]*dXP0mN[IDX_XPs2]
			+ m_np[3]*dXP0mN[IDX_XPs1]
			+ m_np[4]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
DIRK54Solver::dPredStateForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	case 5:
		return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xs4]
			+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A4][DIFFERENTIAL]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A5][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B4][DIFFERENTIAL]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B5][DIFFERENTIAL]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
DIRK54Solver::dPredDerAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return dXP0mN[IDX_XPm1];

	case 2:
		return dXP0mN[IDX_XPs1]; 

	case 3:
		return dXP0mN[IDX_XPs2];

	case 4:
		return dXP0mN[IDX_XPs3];

	case 5:
		return dXP0mN[IDX_XPs4];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal
DIRK54Solver::dPredStateAlgForStageS(unsigned uStage,
	const doublereal dXm1mN[5],
	const doublereal dXP0mN[6]) const
{
	switch (uStage) {
	case 1:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 2:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 3:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 4:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A4][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B4][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	case 5:
		return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xs4]
			+ m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xs3]
			+ m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xs2]
			+ m_a[IDX_A4][ALGEBRAIC]*dXm1mN[IDX_Xs1]
			+ m_a[IDX_A5][ALGEBRAIC]*dXm1mN[IDX_Xm1]
			+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0]
			+ m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPs4]
			+ m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPs3]
			+ m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPs2]
			+ m_b[IDX_B4][ALGEBRAIC]*dXP0mN[IDX_XPs1]
			+ m_b[IDX_B5][ALGEBRAIC]*dXP0mN[IDX_XPm1];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* DIRK54Solver - end */