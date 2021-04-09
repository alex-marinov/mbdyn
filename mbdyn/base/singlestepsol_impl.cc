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

	m_b[IDX_B0][DIFFERENTIAL][IDX_Real] = m_gamma0 * m_gamma2 * dT / m_gamma1;
	m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real] = m_gamma0 * (1. - m_gamma2) * dT / m_gamma1;
	m_b[IDX_BI1][DIFFERENTIAL][IDX_Real] = (m_gamma1 - m_gamma0) * dT / m_gamma1;

	DEBUGCOUT("Predict" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL][IDX_Real] << std::endl
			<< "b1    = " << m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real] << std::endl
			<< "b2    = " << m_b[IDX_BI1][DIFFERENTIAL][IDX_Real] << std::endl);

	m_b[IDX_B0][ALGEBRAIC][IDX_Real] = m_b[IDX_B0][DIFFERENTIAL][IDX_Real];
	m_b[IDX_Bm1][ALGEBRAIC][IDX_Real] = m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real];
	m_b[IDX_BI1][ALGEBRAIC][IDX_Real] = m_b[IDX_BI1][DIFFERENTIAL][IDX_Real];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC][IDX_Real] << std::endl
		<< "b1    = " << m_b[IDX_Bm1][ALGEBRAIC][IDX_Real] << std::endl
		<< "b2    = " << m_b[IDX_BI1][ALGEBRAIC][IDX_Real] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL][IDX_Real];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC][IDX_Real];

	//coefficients for updating intermediate variables
    m_c[IDX_C1][IDX_Real] = - (1. - m_gamma1)/m_gamma1;
    m_d[IDX_D1][IDX_Real] = m_gamma2/m_gamma1;
    m_e[IDX_E1][IDX_Real] = (1. - m_gamma2)/m_gamma1;
}

doublereal
SS2Solver::dPredDer(const doublereal dXm1,
		const doublereal dXP0mI[3][2]) const
{
	return dXP0mI[IDX_XPm1][IDX_Real]; 
}

doublereal
SS2Solver::dPredState(const doublereal dXm1,
		const doublereal dXP0mI[3][2]) const
{
	return dXm1 + m_b[IDX_B0][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XP0][IDX_Real] 
        + m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XPm1][IDX_Real]
        + m_b[IDX_BI1][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XPI1][IDX_Real];
}

doublereal
SS2Solver::dPredDerAlg(const doublereal dXm1,
		const doublereal dXP0mI[3][2]) const
{
	return dXP0mI[IDX_XPm1][IDX_Real]; 
}

doublereal
SS2Solver::dPredStateAlg(const doublereal dXm1,
		const doublereal dXP0mI[3][2]) const
{
	return dXm1 + m_b[IDX_B0][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XP0][IDX_Real]
        + m_b[IDX_Bm1][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XPm1][IDX_Real]
        + m_b[IDX_BI1][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XPI1][IDX_Real];
}

doublereal 
SS2Solver::dUpdateInteReal(unsigned uNumber, const doublereal dXP1[2],
		const doublereal dXP2[2], const doublereal dXP3[2]) const
{
    return m_c[IDX_C1][IDX_Real]*dXP1[IDX_Real] + m_d[IDX_D1][IDX_Real]*dXP2[IDX_Real] + m_e[IDX_E1][IDX_Real]*dXP3[IDX_Real];
}

doublereal 
SS2Solver::dUpdateInteImag(unsigned uNumber, const doublereal dXP1[2],
		const doublereal dXP2[2], const doublereal dXP3[2]) const
{
    return 0.;
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

	doublereal m_gamma0 = 1. / (1. + m_dRho);
	doublereal m_gamma1real = (5. - m_dRho) / (4. * (1. + m_dRho));
	doublereal m_gamma1imag = sqrt(15.) * (1. - m_dRho) / (12. * (1. + m_dRho));
	doublereal m_gamma2 = m_gamma0;
	doublereal m_gamma3real = m_gamma1real;
	doublereal m_gamma3imag = -m_gamma1imag;
	doublereal m_gamma4 = m_gamma2;

	m_b[IDX_B0][DIFFERENTIAL][IDX_Real] = m_gamma0 * m_gamma2 * m_gamma4 * dT / (m_gamma1real * m_gamma3real - m_gamma1imag * m_gamma3imag);
	m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real] = m_gamma0 * m_gamma2 * (1. - m_gamma4) * dT / (m_gamma1real * m_gamma3real - m_gamma1imag * m_gamma3imag);
	m_b[IDX_BI1][DIFFERENTIAL][IDX_Real] = m_gamma0 * (m_gamma3real - m_gamma2) * dT / (m_gamma1real * m_gamma3real - m_gamma1imag * m_gamma3imag);
	m_b[IDX_BI1][DIFFERENTIAL][IDX_Imag] = m_gamma0 * m_gamma3imag * dT / (m_gamma1real * m_gamma3real - m_gamma1imag * m_gamma3imag);
	m_b[IDX_BI2][DIFFERENTIAL][IDX_Real] = (m_gamma1real * m_gamma1real + m_gamma1imag * m_gamma1imag - m_gamma0 * m_gamma1real) * dT / (m_gamma1real * m_gamma1real + m_gamma1imag * m_gamma1imag);
	m_b[IDX_BI2][DIFFERENTIAL][IDX_Imag] = (m_gamma0 * m_gamma1imag) * dT / (m_gamma1real * m_gamma1real + m_gamma1imag * m_gamma1imag);

	DEBUGCOUT("Predict" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL][IDX_Real] << std::endl
			<< "b1    = " << m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real] << std::endl
			<< "b2    = " << m_b[IDX_BI1][DIFFERENTIAL][IDX_Real] << " + " << m_b[IDX_BI1][DIFFERENTIAL][IDX_Imag] << "*i" << std::endl
			<< "b3    = " << m_b[IDX_BI2][DIFFERENTIAL][IDX_Real] << " + " << m_b[IDX_BI2][DIFFERENTIAL][IDX_Imag] << "*i" << std::endl);

	m_b[IDX_B0][ALGEBRAIC][IDX_Real] = m_b[IDX_B0][DIFFERENTIAL][IDX_Real];
	m_b[IDX_Bm1][ALGEBRAIC][IDX_Real] = m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real];
	m_b[IDX_BI1][ALGEBRAIC][IDX_Real] = m_b[IDX_BI1][DIFFERENTIAL][IDX_Real];
	m_b[IDX_BI1][ALGEBRAIC][IDX_Imag] = m_b[IDX_BI1][DIFFERENTIAL][IDX_Imag];
	m_b[IDX_BI2][ALGEBRAIC][IDX_Real] = m_b[IDX_BI2][DIFFERENTIAL][IDX_Real];
	m_b[IDX_BI2][ALGEBRAIC][IDX_Imag] = m_b[IDX_BI2][DIFFERENTIAL][IDX_Imag];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC][IDX_Real] << std::endl
		<< "b1    = " << m_b[IDX_Bm1][ALGEBRAIC][IDX_Real] << std::endl
		<< "b2    = " << m_b[IDX_BI1][ALGEBRAIC][IDX_Real] << " + " << m_b[IDX_BI1][ALGEBRAIC][IDX_Imag] << "*i"  << std::endl
		<< "b3    = " << m_b[IDX_BI2][ALGEBRAIC][IDX_Real] << " + " << m_b[IDX_BI2][ALGEBRAIC][IDX_Imag] << "*i"  << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL][IDX_Real];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC][IDX_Real];

    //coefficients for updating intermediate variables
	m_c[IDX_C1][IDX_Real] = (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag - m_gamma3real) / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);
	m_c[IDX_C1][IDX_Imag] = m_gamma3imag / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);
	m_d[IDX_D1][IDX_Real] = m_gamma4 * m_gamma3real / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);
	m_d[IDX_D1][IDX_Imag] = -m_gamma4 * m_gamma3imag / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);
	m_e[IDX_E1][IDX_Real] = (1. - m_gamma4) * m_gamma3real / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);
	m_e[IDX_E1][IDX_Imag] = -(1. - m_gamma4) * m_gamma3imag / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);	
	m_c[IDX_C2][IDX_Real] = (m_gamma1real * m_gamma1real + m_gamma1imag * m_gamma1imag - m_gamma1real) / (m_gamma1real * m_gamma1real + m_gamma1imag * m_gamma1imag);
	m_c[IDX_C2][IDX_Imag] = m_gamma1imag / (m_gamma1real * m_gamma1real + m_gamma1imag * m_gamma1imag);
	m_d[IDX_D2][IDX_Real] = m_gamma2 * m_gamma1real / (m_gamma1real * m_gamma1real + m_gamma1imag * m_gamma1imag);
	m_d[IDX_D2][IDX_Imag] = -m_gamma2 * m_gamma1imag / (m_gamma1real * m_gamma1real + m_gamma1imag * m_gamma1imag);
	m_e[IDX_E2][IDX_Real] = (1. - m_gamma2) * m_gamma1real / (m_gamma1real * m_gamma1real + m_gamma1imag * m_gamma1imag);
	m_e[IDX_E2][IDX_Imag] = -(1. - m_gamma2) * m_gamma1imag / (m_gamma1real * m_gamma1real + m_gamma1imag * m_gamma1imag);

}

doublereal
SS3Solver::dPredDer(const doublereal dXm1,
		const doublereal dXP0mI[4][2]) const
{
	return dXP0mI[IDX_XPm1][IDX_Real]; 
}

doublereal
SS3Solver::dPredState(const doublereal dXm1,
		const doublereal dXP0mI[4][2]) const
{
	return dXm1 + m_b[IDX_B0][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XP0][IDX_Real] 
        + m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XPm1][IDX_Real]
        + m_b[IDX_BI1][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XPI1][IDX_Real]
		- m_b[IDX_BI1][DIFFERENTIAL][IDX_Imag]*dXP0mI[IDX_XPI1][IDX_Imag]
		+ m_b[IDX_BI2][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XPI2][IDX_Real]
		- m_b[IDX_BI2][DIFFERENTIAL][IDX_Imag]*dXP0mI[IDX_XPI2][IDX_Imag];
}

doublereal
SS3Solver::dPredDerAlg(const doublereal dXm1,
		const doublereal dXP0mI[4][2]) const
{
	return dXP0mI[IDX_XPm1][IDX_Real]; 
}

doublereal
SS3Solver::dPredStateAlg(const doublereal dXm1,
		const doublereal dXP0mI[4][2]) const
{
	return dXm1 + m_b[IDX_B0][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XP0][IDX_Real]
        + m_b[IDX_Bm1][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XPm1][IDX_Real]
        + m_b[IDX_BI1][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XPI1][IDX_Real]
		- m_b[IDX_BI1][ALGEBRAIC][IDX_Imag]*dXP0mI[IDX_XPI1][IDX_Imag]
		+ m_b[IDX_BI2][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XPI2][IDX_Real]
		- m_b[IDX_BI2][ALGEBRAIC][IDX_Imag]*dXP0mI[IDX_XPI2][IDX_Imag];
}

doublereal 
SS3Solver::dUpdateInteReal(unsigned uNumber, const doublereal dXP1[2],
		const doublereal dXP2[2], const doublereal dXP3[2]) const
{
	switch (uNumber)
	{
	case 1:
		return m_c[IDX_C1][IDX_Real]*dXP1[IDX_Real] - m_c[IDX_C1][IDX_Imag]*dXP1[IDX_Imag] 
			+ m_d[IDX_D1][IDX_Real]*dXP2[IDX_Real] - m_d[IDX_D1][IDX_Imag]*dXP2[IDX_Imag]
			+ m_e[IDX_E1][IDX_Real]*dXP3[IDX_Real] - m_e[IDX_E1][IDX_Imag]*dXP3[IDX_Imag];

	case 2:
		return m_c[IDX_C2][IDX_Real]*dXP1[IDX_Real] - m_c[IDX_C2][IDX_Imag]*dXP1[IDX_Imag]
			+ m_d[IDX_D2][IDX_Real]*dXP2[IDX_Real] - m_d[IDX_D2][IDX_Imag]*dXP2[IDX_Imag]  
			+ m_e[IDX_E2][IDX_Real]*dXP3[IDX_Real] - m_e[IDX_E2][IDX_Imag]*dXP3[IDX_Imag];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal 
SS3Solver::dUpdateInteImag(unsigned uNumber, const doublereal dXP1[2],
		const doublereal dXP2[2], const doublereal dXP3[2]) const
{
	switch (uNumber)
	{
	case 1:
		return m_c[IDX_C1][IDX_Real]*dXP1[IDX_Imag] + m_c[IDX_C1][IDX_Imag]*dXP1[IDX_Real] 
			+ m_d[IDX_D1][IDX_Real]*dXP2[IDX_Imag] + m_d[IDX_D1][IDX_Imag]*dXP2[IDX_Real]
			+ m_e[IDX_E1][IDX_Real]*dXP3[IDX_Imag] + m_e[IDX_E1][IDX_Imag]*dXP3[IDX_Real];

	case 2:
		return m_c[IDX_C2][IDX_Real]*dXP1[IDX_Imag] + m_c[IDX_C2][IDX_Imag]*dXP1[IDX_Real]
			+ m_d[IDX_D2][IDX_Real]*dXP2[IDX_Imag] + m_d[IDX_D2][IDX_Imag]*dXP2[IDX_Real]  
			+ m_e[IDX_E2][IDX_Real]*dXP3[IDX_Imag] + m_e[IDX_E2][IDX_Imag]*dXP3[IDX_Real];

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

	doublereal m_gamma0 = 1. / (1. + m_dRho);
	doublereal m_gamma1 = -(4625882646065207. * m_dRho - 18136681528176695.) / (13510798882111488. * (1. + m_dRho));
	doublereal m_gamma2 = m_gamma0;
	doublereal m_gamma3real = -(2129516794990537. * m_dRho - 29151114559213513.) / (27021597764222976. * (1 + m_dRho));
	doublereal m_gamma3imag = -5833664606412816. * sqrt(3.) * (m_dRho - 1.) / (27021597764222976. * (m_dRho + 1.));
	doublereal m_gamma4 = m_gamma2;
	doublereal m_gamma5real = m_gamma3real;
	doublereal m_gamma5imag = -m_gamma3imag;
	doublereal m_gamma6 = m_gamma4;

	m_b[IDX_B0][DIFFERENTIAL][IDX_Real] = m_gamma0 * m_gamma2 * m_gamma4 * m_gamma6 * dT / (m_gamma1 * (m_gamma3real * m_gamma5real - m_gamma3imag * m_gamma5imag));
	m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real] = m_gamma0 * m_gamma2 * m_gamma4 * (1. - m_gamma6) * dT / (m_gamma1 * (m_gamma3real * m_gamma5real - m_gamma3imag * m_gamma5imag));
	m_b[IDX_BI1][DIFFERENTIAL][IDX_Real] = m_gamma0 * m_gamma2 * (m_gamma5real - m_gamma4) * dT / (m_gamma1 * (m_gamma3real * m_gamma5real - m_gamma3imag * m_gamma5imag));
	m_b[IDX_BI1][DIFFERENTIAL][IDX_Imag] = m_gamma0 * m_gamma2 * m_gamma5imag * dT / (m_gamma1 * (m_gamma3real * m_gamma5real - m_gamma3imag * m_gamma5imag));
	m_b[IDX_BI2][DIFFERENTIAL][IDX_Real] = m_gamma0 * (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag - m_gamma2 * m_gamma3real) * dT / (m_gamma1 * (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag));
	m_b[IDX_BI2][DIFFERENTIAL][IDX_Imag] = m_gamma0 * m_gamma2 * m_gamma3imag * dT / (m_gamma1 * (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag));
	m_b[IDX_BI3][DIFFERENTIAL][IDX_Real] = (m_gamma1 - m_gamma0) * dT / m_gamma1;
	m_b[IDX_BI3][DIFFERENTIAL][IDX_Imag] = 0.;

	DEBUGCOUT("Predict" << std::endl
			<< "Alpha = " << dAlpha << std::endl
			<< "Differential coefficients:" << std::endl
			<< "Asymptotic rho =" << m_dRho << std::endl 
			<< "b0    = " << m_b[IDX_B0][DIFFERENTIAL][IDX_Real] << std::endl
			<< "b1    = " << m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real] << std::endl
			<< "b2    = " << m_b[IDX_BI1][DIFFERENTIAL][IDX_Real] << " + " << m_b[IDX_BI1][DIFFERENTIAL][IDX_Imag] << "*i"  << std::endl
			<< "b3    = " << m_b[IDX_BI2][DIFFERENTIAL][IDX_Real] << " + " << m_b[IDX_BI2][DIFFERENTIAL][IDX_Imag] << "*i"  << std::endl
			<< "b4    = " << m_b[IDX_BI3][DIFFERENTIAL][IDX_Real] << std::endl);

	m_b[IDX_B0][ALGEBRAIC][IDX_Real] = m_b[IDX_B0][DIFFERENTIAL][IDX_Real];
	m_b[IDX_Bm1][ALGEBRAIC][IDX_Real] = m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real];
	m_b[IDX_BI1][ALGEBRAIC][IDX_Real] = m_b[IDX_BI1][DIFFERENTIAL][IDX_Real];
	m_b[IDX_BI1][ALGEBRAIC][IDX_Imag] = m_b[IDX_BI1][DIFFERENTIAL][IDX_Imag];
	m_b[IDX_BI2][ALGEBRAIC][IDX_Real] = m_b[IDX_BI2][DIFFERENTIAL][IDX_Real];
	m_b[IDX_BI2][ALGEBRAIC][IDX_Imag] = m_b[IDX_BI2][DIFFERENTIAL][IDX_Imag];
	m_b[IDX_BI3][ALGEBRAIC][IDX_Real] = m_b[IDX_BI3][DIFFERENTIAL][IDX_Real];

	DEBUGCOUT("Algebraic coefficients:" << std::endl
		<< "Asymptotic rho =" << m_dAlgebraicRho << " (ignored)" << std::endl 
		<< "b0    = " << m_b[IDX_B0][ALGEBRAIC][IDX_Real] << std::endl
		<< "b1    = " << m_b[IDX_Bm1][ALGEBRAIC][IDX_Real] << std::endl
		<< "b2    = " << m_b[IDX_BI1][ALGEBRAIC][IDX_Real] << " + " << m_b[IDX_BI1][ALGEBRAIC][IDX_Imag] << "*i"   << std::endl
		<< "b3    = " << m_b[IDX_BI2][ALGEBRAIC][IDX_Real] << " + " << m_b[IDX_BI2][ALGEBRAIC][IDX_Imag] << "*i"   << std::endl
		<< "b4    = " << m_b[IDX_BI3][ALGEBRAIC][IDX_Real] << std::endl);

	db0Differential = m_b[IDX_B0][DIFFERENTIAL][IDX_Real];
	db0Algebraic = m_b[IDX_B0][ALGEBRAIC][IDX_Real];

    //coefficients for updating intermediate variables
	m_c[IDX_C1][IDX_Real] = (m_gamma5real * m_gamma5real + m_gamma5imag * m_gamma5imag - m_gamma5real) / (m_gamma5real * m_gamma5real + m_gamma5imag * m_gamma5imag);
	m_c[IDX_C1][IDX_Imag] = m_gamma5imag / (m_gamma5real * m_gamma5real + m_gamma5imag * m_gamma5imag);
	m_d[IDX_D1][IDX_Real] = m_gamma6 * m_gamma5real / (m_gamma5real * m_gamma5real + m_gamma5imag * m_gamma5imag);
	m_d[IDX_D1][IDX_Imag] = -m_gamma6 * m_gamma5imag / (m_gamma5real * m_gamma5real + m_gamma5imag * m_gamma5imag);
	m_e[IDX_E1][IDX_Real] = (1. - m_gamma6) * m_gamma5real / (m_gamma5real * m_gamma5real + m_gamma5imag * m_gamma5imag);
	m_e[IDX_E1][IDX_Imag] = -(1. - m_gamma6) * m_gamma5imag / (m_gamma5real * m_gamma5real + m_gamma5imag * m_gamma5imag);
	m_c[IDX_C2][IDX_Real] = (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag - m_gamma3real) / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);
	m_c[IDX_C2][IDX_Imag] = m_gamma3imag / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);
	m_d[IDX_D2][IDX_Real] = m_gamma4 * m_gamma3real / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);
	m_d[IDX_D2][IDX_Imag] = -m_gamma4 * m_gamma3imag / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);
	m_e[IDX_E2][IDX_Real] = (1. - m_gamma4) * m_gamma3real / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);
	m_e[IDX_E2][IDX_Imag] = -(1. - m_gamma4) * m_gamma3imag / (m_gamma3real * m_gamma3real + m_gamma3imag * m_gamma3imag);
	m_c[IDX_C3][IDX_Real] = -(1. - m_gamma1) / m_gamma1;
	m_c[IDX_C3][IDX_Imag] = 0.;
	m_d[IDX_D3][IDX_Real] = m_gamma2 / m_gamma1;
	m_d[IDX_D3][IDX_Imag] = 0.;
	m_e[IDX_E3][IDX_Real] = (1. - m_gamma2) / m_gamma1;
	m_e[IDX_E3][IDX_Imag] = 0.;

}

doublereal
SS4Solver::dPredDer(const doublereal dXm1,
		const doublereal dXP0mI[5][2]) const
{
	return dXP0mI[IDX_XPm1][IDX_Real]; 
}

doublereal
SS4Solver::dPredState(const doublereal dXm1,
		const doublereal dXP0mI[5][2]) const
{
	return dXm1 + m_b[IDX_B0][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XP0][IDX_Real] 
        + m_b[IDX_Bm1][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XPm1][IDX_Real]
        + m_b[IDX_BI1][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XPI1][IDX_Real]
		- m_b[IDX_BI1][DIFFERENTIAL][IDX_Imag]*dXP0mI[IDX_XPI1][IDX_Imag]
		+ m_b[IDX_BI2][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XPI2][IDX_Real]
		- m_b[IDX_BI2][DIFFERENTIAL][IDX_Imag]*dXP0mI[IDX_XPI2][IDX_Imag]
		+ m_b[IDX_BI3][DIFFERENTIAL][IDX_Real]*dXP0mI[IDX_XPI3][IDX_Real]
		- m_b[IDX_BI3][DIFFERENTIAL][IDX_Imag]*dXP0mI[IDX_XPI3][IDX_Imag];
}

doublereal
SS4Solver::dPredDerAlg(const doublereal dXm1,
		const doublereal dXP0mI[5][2]) const
{
	return dXP0mI[IDX_XPm1][IDX_Real]; 
}

doublereal
SS4Solver::dPredStateAlg(const doublereal dXm1,
		const doublereal dXP0mI[5][2]) const
{
	return dXm1 + m_b[IDX_B0][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XP0][IDX_Real]
        + m_b[IDX_Bm1][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XPm1][IDX_Real]
        + m_b[IDX_BI1][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XPI1][IDX_Real]
		- m_b[IDX_BI1][ALGEBRAIC][IDX_Imag]*dXP0mI[IDX_XPI1][IDX_Imag]
		+ m_b[IDX_BI2][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XPI2][IDX_Real]
		- m_b[IDX_BI2][ALGEBRAIC][IDX_Imag]*dXP0mI[IDX_XPI2][IDX_Imag]
		+ m_b[IDX_BI3][ALGEBRAIC][IDX_Real]*dXP0mI[IDX_XPI3][IDX_Real]
		- m_b[IDX_BI3][ALGEBRAIC][IDX_Imag]*dXP0mI[IDX_XPI3][IDX_Imag];
}

doublereal 
SS4Solver::dUpdateInteReal(unsigned uNumber, const doublereal dXP1[2],
		const doublereal dXP2[2], const doublereal dXP3[2]) const
{
	switch (uNumber)
	{
	case 1:
		return m_c[IDX_C1][IDX_Real]*dXP1[IDX_Real] - m_c[IDX_C1][IDX_Imag]*dXP1[IDX_Imag]
			+ m_d[IDX_D1][IDX_Real]*dXP2[IDX_Real] - m_d[IDX_D1][IDX_Imag]*dXP2[IDX_Imag]
			+ m_e[IDX_E1][IDX_Real]*dXP3[IDX_Real] - m_e[IDX_E1][IDX_Imag]*dXP3[IDX_Imag];

	case 2:
		return m_c[IDX_C2][IDX_Real]*dXP1[IDX_Real] - m_c[IDX_C2][IDX_Imag]*dXP1[IDX_Imag]
			+ m_d[IDX_D2][IDX_Real]*dXP2[IDX_Real] - m_d[IDX_D2][IDX_Imag]*dXP2[IDX_Imag]
			+ m_e[IDX_E2][IDX_Real]*dXP3[IDX_Real] - m_e[IDX_E2][IDX_Imag]*dXP3[IDX_Imag];

	case 3:
		return m_c[IDX_C3][IDX_Real]*dXP1[IDX_Real] - m_c[IDX_C3][IDX_Imag]*dXP1[IDX_Imag]
			+ m_d[IDX_D3][IDX_Real]*dXP2[IDX_Real] - m_d[IDX_D3][IDX_Imag]*dXP2[IDX_Imag]
			+ m_e[IDX_E3][IDX_Real]*dXP3[IDX_Real] - m_e[IDX_E3][IDX_Imag]*dXP3[IDX_Imag];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

doublereal 
SS4Solver::dUpdateInteImag(unsigned uNumber, const doublereal dXP1[2],
		const doublereal dXP2[2], const doublereal dXP3[2]) const
{
	switch (uNumber)
	{
	case 1:
		return m_c[IDX_C1][IDX_Real]*dXP1[IDX_Imag] + m_c[IDX_C1][IDX_Imag]*dXP1[IDX_Real]
			+ m_d[IDX_D1][IDX_Real]*dXP2[IDX_Imag] + m_d[IDX_D1][IDX_Imag]*dXP2[IDX_Real]
			+ m_e[IDX_E1][IDX_Real]*dXP3[IDX_Imag] + m_e[IDX_E1][IDX_Imag]*dXP3[IDX_Real];

	case 2:
		return m_c[IDX_C2][IDX_Real]*dXP1[IDX_Imag] + m_c[IDX_C2][IDX_Imag]*dXP1[IDX_Real]
			+ m_d[IDX_D2][IDX_Real]*dXP2[IDX_Imag] + m_d[IDX_D2][IDX_Imag]*dXP2[IDX_Real]
			+ m_e[IDX_E2][IDX_Real]*dXP3[IDX_Imag] + m_e[IDX_E2][IDX_Imag]*dXP3[IDX_Real];

	case 3:
		return m_c[IDX_C3][IDX_Real]*dXP1[IDX_Imag] + m_c[IDX_C3][IDX_Imag]*dXP1[IDX_Real]
			+ m_d[IDX_D3][IDX_Real]*dXP2[IDX_Imag] + m_d[IDX_D3][IDX_Imag]*dXP2[IDX_Real]
			+ m_e[IDX_E3][IDX_Real]*dXP3[IDX_Imag] + m_e[IDX_E3][IDX_Imag]*dXP3[IDX_Real];

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* SS4Solver - end */
