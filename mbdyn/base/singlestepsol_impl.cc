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
#include "singlestepsol_impl.h"
#include "stepsol.hc"

/* SS2Solver - begin */

SS2Solver::SS2Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	const DriveCaller* pRho,
	const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplSingleStepIntegrator<2>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

SS2Solver::~SS2Solver(void)
{
	NO_OP;
}

void
SS2Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
SS2Solver::SetCoef(doublereal dT, 
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
    m_dRho = m_Rho.dGet();
    m_dAlgebraicRho = m_AlgebraicRho.dGet();

    doublereal m_gamma0 = 1./(1. + m_dRho);
    doublereal m_gamma1 = (3. - m_dRho)/(2.*(1. + m_dRho));
    doublereal m_gamma2 = m_gamma0;

    m_b[IDX_B0][DIFFERENTIAL] = m_gamma0*m_gamma2*dT/m_gamma1;
    m_b[IDX_Bm1][DIFFERENTIAL] = m_gamma0*(1. - m_gamma2)*dT/m_gamma1;
    m_b[IDX_BI1][DIFFERENTIAL] = (m_gamma1 - m_gamma0)*dT/m_gamma1;

    DEBUGCOUT("Predict" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_Bm1][DIFFERENTIAL] << std::endl
			<< "b2    = " << m_b[IDX_BI1][DIFFERENTIAL] << std::endl);

	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_Bm1][ALGEBRAIC] = m_b[IDX_Bm1][DIFFERENTIAL];
	m_b[IDX_BI1][ALGEBRAIC] = m_b[IDX_BI1][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_Bm1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_BI1][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];

    //coefficients for updating intermediate variables
    m_c[IDX_C1] = - (1. - m_gamma1)/m_gamma1;
    m_d[IDX_D1] = m_gamma2/m_gamma1;
    m_e[IDX_E1] = (1. - m_gamma2)/m_gamma1;
}

doublereal
SS2Solver::dPredDer(const doublereal dXm1,
		const doublereal dXP0mI[3]) const
{
	return dXP0mI[IDX_XPm1]; 
}

doublereal
SS2Solver::dPredState(const doublereal dXm1,
		const doublereal dXP0mI[3]) const
{
	return dXm1 + m_b[IDX_B0][DIFFERENTIAL]*dXP0mI[IDX_XP0] 
        + m_b[IDX_Bm1][DIFFERENTIAL]*dXP0mI[IDX_XPm1]
        + m_b[IDX_BI1][DIFFERENTIAL]*dXP0mI[IDX_XPI1];
}

doublereal
SS2Solver::dPredDerAlg(const doublereal dXm1,
		const doublereal dXP0mI[3]) const
{
	return dXP0mI[IDX_XPm1]; 
}

doublereal
SS2Solver::dPredStateAlg(const doublereal dXm1,
		const doublereal dXP0mI[3]) const
{
	return dXm1 + m_b[IDX_B0][ALGEBRAIC]*dXP0mI[IDX_XP0]
        + m_b[IDX_Bm1][ALGEBRAIC]*dXP0mI[IDX_XPm1]
        + m_b[IDX_BI1][ALGEBRAIC]*dXP0mI[IDX_XPI1];
}

doublereal 
SS2Solver::dUpdateInte(unsigned uNumber, const doublereal dXPmI,
		const doublereal dXPIm, const doublereal dXPmIm) const
{
    return m_c[IDX_C1]*dXPmI + m_d[IDX_D1]*dXPIm + m_e[IDX_E1]*dXPmIm;
}

/* SS2Solver - end */

/* SS3Solver - begin */

SS3Solver::SS3Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	const DriveCaller* pRho,
	const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplSingleStepIntegrator<3>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

SS3Solver::~SS3Solver(void)
{
	NO_OP;
}

void
SS3Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
SS3Solver::SetCoef(doublereal dT, 
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
    m_dRho = m_Rho.dGet();
    m_dAlgebraicRho = m_AlgebraicRho.dGet();

    doublereal m_gamma0 = 1./(1. + m_dRho);
    doublereal m_gamma1 = (5. - m_dRho)/(4.*(1. + m_dRho));
    doublereal m_gamma2 = m_gamma0;
	doublereal m_gamma3 = m_gamma1;
	doublereal m_gamma4 = m_gamma2;

    m_b[IDX_B0][DIFFERENTIAL] = m_gamma0*m_gamma2*m_gamma4*dT/(m_gamma1*m_gamma3);
    m_b[IDX_Bm1][DIFFERENTIAL] = m_gamma0*m_gamma2*(1. - m_gamma4)*dT/(m_gamma1*m_gamma3);
    m_b[IDX_BI1][DIFFERENTIAL] = m_gamma0*(m_gamma3 - m_gamma2)*dT/(m_gamma1*m_gamma3);
	m_b[IDX_BI2][DIFFERENTIAL] = (m_gamma1 - m_gamma0)*dT/m_gamma1;

    DEBUGCOUT("Predict" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_Bm1][DIFFERENTIAL] << std::endl
			<< "b2    = " << m_b[IDX_BI1][DIFFERENTIAL] << std::endl
			<< "b3    = " << m_b[IDX_BI2][DIFFERENTIAL] << std::endl);

	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_Bm1][ALGEBRAIC] = m_b[IDX_Bm1][DIFFERENTIAL];
	m_b[IDX_BI1][ALGEBRAIC] = m_b[IDX_BI1][DIFFERENTIAL];
	m_b[IDX_BI2][ALGEBRAIC] = m_b[IDX_BI2][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_Bm1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_BI1][ALGEBRAIC] << std::endl
		<< "b3    = " << m_b[IDX_BI2][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];

    //coefficients for updating intermediate variables
	m_c[IDX_C1] = - (1. - m_gamma3)/m_gamma3;
    m_d[IDX_D1] = m_gamma4/m_gamma3;
    m_e[IDX_E1] = (1. - m_gamma4)/m_gamma3;
    m_c[IDX_C2] = - (1. - m_gamma1)/m_gamma1;
    m_d[IDX_D2] = m_gamma2/m_gamma1;
    m_e[IDX_E2] = (1. - m_gamma2)/m_gamma1;
}

doublereal
SS3Solver::dPredDer(const doublereal dXm1,
		const doublereal dXP0mI[4]) const
{
	return dXP0mI[IDX_XPm1]; 
}

doublereal
SS3Solver::dPredState(const doublereal dXm1,
		const doublereal dXP0mI[4]) const
{
	return dXm1 + m_b[IDX_B0][DIFFERENTIAL]*dXP0mI[IDX_XP0] 
        + m_b[IDX_Bm1][DIFFERENTIAL]*dXP0mI[IDX_XPm1]
        + m_b[IDX_BI1][DIFFERENTIAL]*dXP0mI[IDX_XPI1]
		+ m_b[IDX_BI2][DIFFERENTIAL]*dXP0mI[IDX_XPI2];
}

doublereal
SS3Solver::dPredDerAlg(const doublereal dXm1,
		const doublereal dXP0mI[4]) const
{
	return dXP0mI[IDX_XPm1]; 
}

doublereal
SS3Solver::dPredStateAlg(const doublereal dXm1,
		const doublereal dXP0mI[4]) const
{
	return dXm1 + m_b[IDX_B0][ALGEBRAIC]*dXP0mI[IDX_XP0]
        + m_b[IDX_Bm1][ALGEBRAIC]*dXP0mI[IDX_XPm1]
        + m_b[IDX_BI1][ALGEBRAIC]*dXP0mI[IDX_XPI1]
		+ m_b[IDX_BI2][ALGEBRAIC]*dXP0mI[IDX_XPI2];
}

doublereal 
SS3Solver::dUpdateInte(unsigned uNumber, const doublereal dXPmI,
		const doublereal dXPIm, const doublereal dXPmIm) const
{
	switch (uNumber)
	{
	case 1:
		return m_c[IDX_C1]*dXPmI + m_d[IDX_D1]*dXPIm + m_e[IDX_E1]*dXPmIm;

	case 2:
		return m_c[IDX_C2]*dXPmI + m_d[IDX_D2]*dXPIm + m_e[IDX_E2]*dXPmIm;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* SS3Solver - end */

/* SS4Solver - begin */

SS4Solver::SS4Solver(const doublereal Tl,
	const doublereal dSolTl,
	const integer iMaxIt,
	const DriveCaller* pRho,
	const DriveCaller* pAlgRho,
	const bool bmod_res_test)
: tplSingleStepIntegrator<4>(iMaxIt, Tl, dSolTl, bmod_res_test),
m_Rho(pRho), m_AlgebraicRho(pAlgRho)
{
	ASSERT(pRho != 0);
	ASSERT(pAlgRho != 0);
}

SS4Solver::~SS4Solver(void)
{
	NO_OP;
}

void
SS4Solver::SetDriveHandler(const DriveHandler* pDH)
{
	m_Rho.pGetDriveCaller()->SetDrvHdl(pDH);
	m_AlgebraicRho.pGetDriveCaller()->SetDrvHdl(pDH);
}

void
SS4Solver::SetCoef(doublereal dT, 
		doublereal dAlpha,
		enum StepChange /* NewStep */)
{
    m_dRho = m_Rho.dGet();
    m_dAlgebraicRho = m_AlgebraicRho.dGet();

    doublereal m_gamma0 = 1./(1. + m_dRho);
    doublereal m_gamma1 = (7. - m_dRho)/(6.*(1. + m_dRho));
    doublereal m_gamma2 = m_gamma0;
	doublereal m_gamma3 = m_gamma1;
	doublereal m_gamma4 = m_gamma2;
	doublereal m_gamma5 = m_gamma3;
	doublereal m_gamma6 = m_gamma4;

    m_b[IDX_B0][DIFFERENTIAL] = m_gamma0*m_gamma2*m_gamma4*m_gamma6*dT/(m_gamma1*m_gamma3*m_gamma5);
    m_b[IDX_Bm1][DIFFERENTIAL] = m_gamma0*m_gamma2*m_gamma4*(1. - m_gamma6)*dT/(m_gamma1*m_gamma3*m_gamma5);
    m_b[IDX_BI1][DIFFERENTIAL] = m_gamma0*m_gamma2*(m_gamma5 - m_gamma4)*dT/(m_gamma1*m_gamma3*m_gamma5);
	m_b[IDX_BI2][DIFFERENTIAL] = m_gamma0*(m_gamma3 - m_gamma2)*dT/(m_gamma1*m_gamma3);
	m_b[IDX_BI3][DIFFERENTIAL] = (m_gamma1 - m_gamma0)*dT/m_gamma1;

    DEBUGCOUT("Predict" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL] << std::endl
			<< "b1    = " << m_b[IDX_Bm1][DIFFERENTIAL] << std::endl
			<< "b2    = " << m_b[IDX_BI1][DIFFERENTIAL] << std::endl
			<< "b3    = " << m_b[IDX_BI2][DIFFERENTIAL] << std::endl
			<< "b4    = " << m_b[IDX_BI3][DIFFERENTIAL] << std::endl);

	m_b[IDX_B0][ALGEBRAIC] = m_b[IDX_B0][DIFFERENTIAL];
	m_b[IDX_Bm1][ALGEBRAIC] = m_b[IDX_Bm1][DIFFERENTIAL];
	m_b[IDX_BI1][ALGEBRAIC] = m_b[IDX_BI1][DIFFERENTIAL];
	m_b[IDX_BI2][ALGEBRAIC] = m_b[IDX_BI2][DIFFERENTIAL];
	m_b[IDX_BI3][ALGEBRAIC] = m_b[IDX_BI3][DIFFERENTIAL];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC] << std::endl
		<< "b1    = " << m_b[IDX_Bm1][ALGEBRAIC] << std::endl
		<< "b2    = " << m_b[IDX_BI1][ALGEBRAIC] << std::endl
		<< "b3    = " << m_b[IDX_BI2][ALGEBRAIC] << std::endl
		<< "b4    = " << m_b[IDX_BI3][ALGEBRAIC] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC];

    //coefficients for updating intermediate variables
	m_c[IDX_C1] = - (1. - m_gamma5)/m_gamma5;
    m_d[IDX_D1] = m_gamma6/m_gamma5;
    m_e[IDX_E1] = (1. - m_gamma6)/m_gamma5;
	m_c[IDX_C2] = - (1. - m_gamma3)/m_gamma3;
    m_d[IDX_D2] = m_gamma4/m_gamma3;
    m_e[IDX_E2] = (1. - m_gamma4)/m_gamma3;
    m_c[IDX_C3] = - (1. - m_gamma1)/m_gamma1;
    m_d[IDX_D3] = m_gamma2/m_gamma1;
    m_e[IDX_E3] = (1. - m_gamma2)/m_gamma1;
}

doublereal
SS4Solver::dPredDer(const doublereal dXm1,
		const doublereal dXP0mI[5]) const
{
	return dXP0mI[IDX_XPm1]; 
}

doublereal
SS4Solver::dPredState(const doublereal dXm1,
		const doublereal dXP0mI[5]) const
{
	return dXm1 + m_b[IDX_B0][DIFFERENTIAL]*dXP0mI[IDX_XP0] 
        + m_b[IDX_Bm1][DIFFERENTIAL]*dXP0mI[IDX_XPm1]
        + m_b[IDX_BI1][DIFFERENTIAL]*dXP0mI[IDX_XPI1]
		+ m_b[IDX_BI2][DIFFERENTIAL]*dXP0mI[IDX_XPI2]
		+ m_b[IDX_BI3][DIFFERENTIAL]*dXP0mI[IDX_XPI3];
}

doublereal
SS4Solver::dPredDerAlg(const doublereal dXm1,
		const doublereal dXP0mI[5]) const
{
	return dXP0mI[IDX_XPm1]; 
}

doublereal
SS4Solver::dPredStateAlg(const doublereal dXm1,
		const doublereal dXP0mI[5]) const
{
	return dXm1 + m_b[IDX_B0][ALGEBRAIC]*dXP0mI[IDX_XP0]
        + m_b[IDX_Bm1][ALGEBRAIC]*dXP0mI[IDX_XPm1]
        + m_b[IDX_BI1][ALGEBRAIC]*dXP0mI[IDX_XPI1]
		+ m_b[IDX_BI2][ALGEBRAIC]*dXP0mI[IDX_XPI2]
		+ m_b[IDX_BI3][ALGEBRAIC]*dXP0mI[IDX_XPI3];
}

doublereal 
SS4Solver::dUpdateInte(unsigned uNumber, const doublereal dXPmI,
		const doublereal dXPIm, const doublereal dXPmIm) const
{
	switch (uNumber)
	{
	case 1:
		return m_c[IDX_C1]*dXPmI + m_d[IDX_D1]*dXPIm + m_e[IDX_E1]*dXPmIm;

	case 2:
		return m_c[IDX_C2]*dXPmI + m_d[IDX_D2]*dXPIm + m_e[IDX_E2]*dXPmIm;

	case 3:
		return m_c[IDX_C3]*dXPmI + m_d[IDX_D3]*dXPIm + m_e[IDX_E3]*dXPmIm;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* SS4Solver - end */