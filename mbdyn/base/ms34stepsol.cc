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
#include "ms34stepsol.h"
#include "stepsol.hc"

/* TunableStep3Solver - begin */

TunableStep3Solver::TunableStep3Solver(const doublereal Tl,
		const doublereal dSolTl,
		const integer iMaxIt,
		const DriveCaller* pRho,
		const DriveCaller* pAlgRho,
		const bool bmod_res_test)
: tplStepNIntegrator<3>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

TunableStep3Solver::~TunableStep3Solver(void)
{
	NO_OP;
}

void
TunableStep3Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
TunableStep3Solver::SetCoef(doublereal dT,
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
	doublereal dRho = m_Rho.dGet();
	doublereal dAlgebraicRho = m_AlgebraicRho.dGet();

	m_mp[0] = -6.*dAlpha*dAlpha*(1. + dAlpha)/dT;
	m_mp[1] = -m_mp[0];
	m_np[0] = 1. + 4.*dAlpha + 3.*dAlpha*dAlpha;
	m_np[1] = dAlpha*(2. + 3.*dAlpha);
	
	/*doublereal dDen = 2.*(1.+dAlpha)-(1.-dRho)*(1.-dRho);
	doublereal dBeta = dAlpha*((1.-dRho)*(1.-dRho)*(2.+dAlpha)
		+2.*(2.*dRho-1.)*(1.+dAlpha))/dDen;
	doublereal dDelta = .5*dAlpha*dAlpha*(1.-dRho)*(1.-dRho)/dDen;


	a[0][DIFFERENTIAL] = 1.-dBeta;
	a[1][DIFFERENTIAL] = dBeta;
	a[2][DIFFERENTIAL] = 0.;
	b[0][DIFFERENTIAL] = dT*(dDelta/dAlpha+dAlpha/2);
	b[1][DIFFERENTIAL] = dT*(dBeta/2.+dAlpha/2.-dDelta/dAlpha*(1.+dAlpha));
	b[2][DIFFERENTIAL] = dT*(dBeta/2.+dDelta);
	b[3][DIFFERENTIAL] = 0.;*/

	doublereal dDen = dRho*dRho - 5.*dRho + 10.;
	doublereal dBeta = (1. + dRho)*dDen;
	
	m_a[IDX_A1][DIFFERENTIAL] = 3.*(2.*dRho*dRho - 9.*dRho + 5.)/dDen;
	m_a[IDX_A2][DIFFERENTIAL] = -3.*(5.*dRho*dRho - 9.*dRho + 2.)/dDen;
	m_a[IDX_A3][DIFFERENTIAL] = (10.*dRho*dRho - 5.*dRho + 1.)/dDen;
	m_b[IDX_B0][DIFFERENTIAL] = dT*(6./dBeta);
	m_b[IDX_B1][DIFFERENTIAL] = dT*(18.*dRho/dBeta);
	m_b[IDX_B2][DIFFERENTIAL] = dT*(18.*dRho*dRho/dBeta);
	m_b[IDX_B3][DIFFERENTIAL] = dT*(6.*dRho*dRho*dRho/dBeta);

	DEBUGCOUT("Predict()" << std::endl
		<< "Alpha = " << dAlpha << std::endl
		<< "Differential coefficients:" << std::endl
		<< "Asymptotic rho =" << dRho << std::endl 
		<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
		<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
		<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
		<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
		<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
		<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
		<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl);

	if (dAlgebraicRho != dRho) {

		/*doublereal dDen = 2.*(1.+dAlpha)-(1.-dAlgebraicRho)*(1.-dAlgebraicRho);
		doublereal dBeta = dAlpha*((1.-dAlgebraicRho)*(1.-dAlgebraicRho)*(2.+dAlpha)
				+2.*(2.*dAlgebraicRho-1.)*(1.+dAlpha))/dDen;
		doublereal dDelta = .5*dAlpha*dAlpha*(1.-dAlgebraicRho)*(1.-dAlgebraicRho)/dDen;

		a[1][ALGEBRAIC] = dBeta;
		b[0][ALGEBRAIC] = dT*(dDelta/dAlpha+dAlpha/2.);
		b[1][ALGEBRAIC] = dT*(dBeta/2.+dAlpha/2.-dDelta/dAlpha*(1.+dAlpha));
		b[2][ALGEBRAIC] = dT*(dBeta/2.+dDelta);*/
		m_a[IDX_A1][ALGEBRAIC] = 3.*(2.*dAlgebraicRho*dAlgebraicRho - 9.*dAlgebraicRho + 5.)/(dAlgebraicRho*dAlgebraicRho - 5.*dAlgebraicRho + 10.);
		m_a[IDX_A2][ALGEBRAIC] = -3.*(5.*dAlgebraicRho*dAlgebraicRho - 9.*dAlgebraicRho + 2.)/(dAlgebraicRho*dAlgebraicRho - 5.*dAlgebraicRho + 10.);
		m_a[IDX_A3][ALGEBRAIC] = (10.*dAlgebraicRho*dAlgebraicRho - 5.*dAlgebraicRho+1.)/(dAlgebraicRho*dAlgebraicRho - 5.*dAlgebraicRho + 10.);
		m_b[IDX_B0][ALGEBRAIC] = dT*(6./((1. + dAlgebraicRho)*(dAlgebraicRho*dAlgebraicRho - 5.*dAlgebraicRho + 10.)));
		m_b[IDX_B1][ALGEBRAIC] = dT*(18.*dAlgebraicRho/((1. + dAlgebraicRho)*(dAlgebraicRho*dAlgebraicRho - 5.*dAlgebraicRho + 10.)));
		m_b[IDX_B2][ALGEBRAIC] = dT*(18.*dAlgebraicRho*dAlgebraicRho/((1. + dAlgebraicRho)*(dAlgebraicRho*dAlgebraicRho - 5.*dAlgebraicRho + 10.)));
		m_b[IDX_B3][ALGEBRAIC] = dT*(6.*dAlgebraicRho*dAlgebraicRho*dAlgebraicRho/((1. + dAlgebraicRho)*(dAlgebraicRho*dAlgebraicRho - 5.*dAlgebraicRho + 10.)));

	} else {

		/*doublereal dDen = 2.*(1.+dAlpha)-(1.-dAlgebraicRho)*(1.-dAlgebraicRho);
		doublereal dBeta = dAlpha*((1.-dAlgebraicRho)*(1.-dAlgebraicRho)*(2.+dAlpha)
				+2.*(2.*dAlgebraicRho-1.)*(1.+dAlpha))/dDen;
		doublereal dDelta = .5*dAlpha*dAlpha*(1.-dAlgebraicRho)*(1.-dAlgebraicRho)/dDen;

		a[1][ALGEBRAIC] = dBeta;
		b[0][ALGEBRAIC] = dT*(dDelta/dAlpha+dAlpha/2.);
		b[1][ALGEBRAIC] = dT*(dBeta/2.+dAlpha/2.-dDelta/dAlpha*(1.+dAlpha));
		b[2][ALGEBRAIC] = dT*(dBeta/2.+dDelta);*/
		m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
		m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
		m_a[IDX_A3][ALGEBRAIC] = m_a[IDX_A3][DIFFERENTIAL];
		m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
		m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
		m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
		m_b[IDX_B3][ALGEBRAIC] = m_b[IDX_B3][DIFFERENTIAL];
	}

	DEBUGCOUT("Algebraic coefficients:" << std::endl
			<< "Asymptotic rho =" << dAlgebraicRho << std::endl 
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
TunableStep3Solver::dPredDer(const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const
{
	/*return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;*/
	return dXP0mN[IDX_XPm1];
}

doublereal
TunableStep3Solver::dPredState(const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const
{
	return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1] + m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm2] + m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xm3]
		+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0] + m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1] + m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm2] + m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPm3];
}

doublereal
TunableStep3Solver::dPredDerAlg(const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const
{
	/*return np[0]*dXPm1+np[1]*dXPm2-mp[1]*dXm1;*/
	return dXP0mN[IDX_XPm1];
}

doublereal
TunableStep3Solver::dPredStateAlg(const doublereal dXm1mN[3],
		const doublereal dXP0mN[4]) const
{
	return m_a[IDX_A1][ALGEBRAIC]*dXm1mN[IDX_Xm1] + m_a[IDX_A2][ALGEBRAIC]*dXm1mN[IDX_Xm2] + m_a[IDX_A3][ALGEBRAIC]*dXm1mN[IDX_Xm3]
		+ m_b[IDX_B0][ALGEBRAIC]*dXP0mN[IDX_XP0] + m_b[IDX_B1][ALGEBRAIC]*dXP0mN[IDX_XPm1] + m_b[IDX_B2][ALGEBRAIC]*dXP0mN[IDX_XPm2] + m_b[IDX_B3][ALGEBRAIC]*dXP0mN[IDX_XPm3];
	/*return dXm1;*/
}

/* TunableStep3Solver - end */


/* Step4Solver - begin */

TunableStep4Solver::TunableStep4Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	const DriveCaller* pRho,
	const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplStepNIntegrator<4>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

TunableStep4Solver::~TunableStep4Solver(void)
{
	NO_OP;
}

void
TunableStep4Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
TunableStep4Solver::SetCoef(doublereal dT,
	doublereal dAlpha,
	enum StepChange /* NewStep */)
{
	doublereal dRho = m_Rho.dGet();
	doublereal dAlgebraicRho = m_AlgebraicRho.dGet();


	/*mp[0] = -6.*dAlpha*dAlpha*(1.+dAlpha)/dT;
	mp[1] = -mp[0];
	np[0] = 1.+4.*dAlpha+3.*dAlpha*dAlpha;
	np[1] = dAlpha*(2.+3.*dAlpha);*/

	m_a[IDX_A1][DIFFERENTIAL] = 4.0*(-2.0*dRho*dRho*dRho + 13.0*dRho*dRho - 35.0*dRho + 14.0)/(-dRho*dRho*dRho + 7.0*dRho*dRho - 21.0*dRho + 35.0);
	m_a[IDX_A2][DIFFERENTIAL] = -4.0*(1.0 - dRho)*(7.0*dRho*dRho - 34.0*dRho + 7.0)/(-dRho*dRho*dRho + 7.0*dRho*dRho - 21.0*dRho + 35.0);
	m_a[IDX_A3][DIFFERENTIAL] = -4.0*(14.0*dRho*dRho*dRho - 35.0*dRho*dRho + 13.0*dRho - 2.0)/(-dRho*dRho*dRho + 7.0*dRho*dRho - 21.0*dRho + 35.0);
	m_a[IDX_A4][DIFFERENTIAL] = (35.0*dRho*dRho*dRho - 21.0*dRho*dRho + 7.0*dRho - 1.0)/(-dRho*dRho*dRho + 7.0*dRho*dRho - 21.0*dRho + 35.0);
	m_b[IDX_B0][DIFFERENTIAL] = dT*(20.0/((1.0 + dRho)*(-dRho*dRho*dRho + 7.0*dRho*dRho - 21.0*dRho + 35.0)));
	m_b[IDX_B1][DIFFERENTIAL] = dT*(80.0*dRho/((1.0 + dRho)*(-dRho*dRho*dRho + 7.0*dRho*dRho - 21.0*dRho + 35.0)));
	m_b[IDX_B2][DIFFERENTIAL] = dT*(120.0*dRho*dRho/((1.0 + dRho)*(-dRho*dRho*dRho + 7.0*dRho*dRho - 21.0*dRho + 35.0)));
	m_b[IDX_B3][DIFFERENTIAL] = dT*(80.0*dRho*dRho*dRho/((1.0 + dRho)*(-dRho*dRho*dRho + 7.0*dRho*dRho - 21.0*dRho + 35.0)));
	m_b[IDX_B4][DIFFERENTIAL] = dT*(20.0*dRho*dRho*dRho*dRho/((1.0 + dRho)*(-dRho*dRho*dRho + 7.0*dRho*dRho - 21.0*dRho + 35.0)));

	DEBUGCOUT("Predict()" << std::endl
		<< "Alpha = " << dAlpha << std::endl
		<< "Differential coefficients:" << std::endl
		<< "Asymptotic rho =" << dRho << std::endl 
		<< "a1    = " << m_a[IDX_A1][DIFFERENTIAL] << std::endl
		<< "a2    = " << m_a[IDX_A2][DIFFERENTIAL] << std::endl
		<< "a3    = " << m_a[IDX_A3][DIFFERENTIAL] << std::endl
		<< "a4    = " << m_a[IDX_A4][DIFFERENTIAL] << std::endl
		<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
		<< "b1    = " << m_b[IDX_B1][DIFFERENTIAL] << std::endl
		<< "b2    = " << m_b[IDX_B2][DIFFERENTIAL] << std::endl
		<< "b3    = " << m_b[IDX_B3][DIFFERENTIAL] << std::endl
		<< "b4    = " << m_b[IDX_B4][DIFFERENTIAL] << std::endl);

	if (dAlgebraicRho != dRho) {
		m_a[IDX_A1][ALGEBRAIC] = 4.0*(-2.0*dAlgebraicRho*dAlgebraicRho*dAlgebraicRho+13.0*dAlgebraicRho*dAlgebraicRho-35.0*dAlgebraicRho+14.0)/(-dAlgebraicRho*dAlgebraicRho*dAlgebraicRho+7.0*dAlgebraicRho*dAlgebraicRho-21.0*dAlgebraicRho+35.0);
		m_a[IDX_A2][ALGEBRAIC] = -4.0*(1.0-dAlgebraicRho)*(7.0*dAlgebraicRho*dAlgebraicRho-34.0*dAlgebraicRho+7.0)/(-dAlgebraicRho*dAlgebraicRho*dAlgebraicRho+7.0*dAlgebraicRho*dAlgebraicRho-21.0*dAlgebraicRho+35.0);
		m_a[IDX_A3][ALGEBRAIC] = -4.0*(14.0*dAlgebraicRho*dAlgebraicRho*dAlgebraicRho-35.0*dAlgebraicRho*dAlgebraicRho+13.0*dAlgebraicRho-2.0)/(-dAlgebraicRho*dAlgebraicRho*dAlgebraicRho+7.0*dAlgebraicRho*dAlgebraicRho-21.0*dAlgebraicRho+35.0);
		m_a[IDX_A4][ALGEBRAIC] = (35.0*dAlgebraicRho*dAlgebraicRho*dAlgebraicRho-21.0*dAlgebraicRho*dAlgebraicRho+7.0*dAlgebraicRho-1.0)/(-dAlgebraicRho*dAlgebraicRho*dAlgebraicRho+7.0*dAlgebraicRho*dAlgebraicRho-21.0*dAlgebraicRho+35.0);
		m_b[IDX_B0][ALGEBRAIC] = dT*(20.0/((1.0+dAlgebraicRho)*(-dAlgebraicRho*dAlgebraicRho*dAlgebraicRho+7.0*dAlgebraicRho*dAlgebraicRho-21.0*dAlgebraicRho+35.0)));
		m_b[IDX_B1][ALGEBRAIC] = dT*(80.0*dAlgebraicRho/((1.0+dAlgebraicRho)*(-dAlgebraicRho*dAlgebraicRho*dAlgebraicRho+7.0*dAlgebraicRho*dAlgebraicRho-21.0*dAlgebraicRho+35.0)));
		m_b[IDX_B2][ALGEBRAIC] = dT*(120.0*dAlgebraicRho*dAlgebraicRho/((1.0+dAlgebraicRho)*(-dAlgebraicRho*dAlgebraicRho*dAlgebraicRho+7.0*dAlgebraicRho*dAlgebraicRho-21.0*dAlgebraicRho+35.0)));
		m_b[IDX_B3][ALGEBRAIC] = dT*(80.0*dAlgebraicRho*dAlgebraicRho*dAlgebraicRho/((1.0+dAlgebraicRho)*(-dAlgebraicRho*dAlgebraicRho*dAlgebraicRho+7.0*dAlgebraicRho*dAlgebraicRho-21.0*dAlgebraicRho+35.0)));
		m_b[IDX_B4][ALGEBRAIC] = dT*(20.0*dAlgebraicRho*dAlgebraicRho*dAlgebraicRho*dAlgebraicRho/((1.0+dAlgebraicRho)*(-dAlgebraicRho*dAlgebraicRho*dAlgebraicRho+7.0*dAlgebraicRho*dAlgebraicRho-21.0*dAlgebraicRho+35.0)));

	} else {
		m_a[IDX_A1][ALGEBRAIC] = m_a[IDX_A1][DIFFERENTIAL];
		m_a[IDX_A2][ALGEBRAIC] = m_a[IDX_A2][DIFFERENTIAL];
		m_a[IDX_A3][ALGEBRAIC] = m_a[IDX_A3][DIFFERENTIAL];
		m_a[IDX_A4][ALGEBRAIC] = m_a[IDX_A4][DIFFERENTIAL];
		m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
		m_b[IDX_B1][ALGEBRAIC] = m_b[IDX_B1][DIFFERENTIAL];
		m_b[IDX_B2][ALGEBRAIC] = m_b[IDX_B2][DIFFERENTIAL];
		m_b[IDX_B3][ALGEBRAIC] = m_b[IDX_B3][DIFFERENTIAL];
		m_b[IDX_B4][ALGEBRAIC] = m_b[IDX_B4][DIFFERENTIAL];
	}

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << dAlgebraicRho << std::endl 
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
TunableStep4Solver::dPredDer(const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const
{
	/*return mp[0]*dXm1+mp[1]*dXm2+np[0]*dXPm1+np[1]*dXPm2;*/
	return dXP0mN[IDX_XPm1];
}

doublereal
TunableStep4Solver::dPredState(const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const
{
	return m_a[IDX_A1][DIFFERENTIAL]*dXm1mN[IDX_Xm1]
		+ m_a[IDX_A2][DIFFERENTIAL]*dXm1mN[IDX_Xm2]
		+ m_a[IDX_A3][DIFFERENTIAL]*dXm1mN[IDX_Xm3]
		+ m_a[IDX_A4][DIFFERENTIAL]*dXm1mN[IDX_Xm4]
		+ m_b[IDX_B0][DIFFERENTIAL]*dXP0mN[IDX_XP0]
		+ m_b[IDX_B1][DIFFERENTIAL]*dXP0mN[IDX_XPm1]
		+ m_b[IDX_B2][DIFFERENTIAL]*dXP0mN[IDX_XPm2]
		+ m_b[IDX_B3][DIFFERENTIAL]*dXP0mN[IDX_XPm3]
		+ m_b[IDX_B4][DIFFERENTIAL]*dXP0mN[IDX_XPm4];
}

doublereal
TunableStep4Solver::dPredDerAlg(const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const
{
	/*return np[0]*dXPm1+np[1]*dXPm2-mp[1]*dXm1;*/
	return dXP0mN[IDX_XPm1];
}

doublereal
TunableStep4Solver::dPredStateAlg(const doublereal dXm1mN[4],
		const doublereal dXP0mN[5]) const
{
	/*return b[0][ALGEBRAIC]*dXP+b[1][ALGEBRAIC]*dXPm1+b[2][ALGEBRAIC]*dXPm2
		-a[1][ALGEBRAIC]*dXm1;*/
	return dXm1mN[IDX_Xm1];
	// FIXME?!?
}

/* Step4Solver - end */
