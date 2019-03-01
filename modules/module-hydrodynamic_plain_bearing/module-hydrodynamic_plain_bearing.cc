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
 AUTHOR: Reinhard Resch <R.RESCH@secop.com>
  Copyright (C) 2011(-2017) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include <limits>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cassert>
#include <cmath>

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>
#include <userelem.h>

#include "module-hydrodynamic_plain_bearing.h"
#include "hydrodynamic_plain_bearing_force.h"

class HydrodynamicPlainBearing: virtual public Elem, public UserDefinedElem
{
public:
    HydrodynamicPlainBearing(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
    virtual ~HydrodynamicPlainBearing(void);
	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		   doublereal dCoef,
		   const VectorHandler& XCurr,
		   const VectorHandler& XPrimeCurr);
	inline
	void ComputeResidual(Vec3& e_R2,
						 Vec3& e_dot_R2,
						 doublereal omega_proj[2],
						 Vec3& F2_R2,
						 Vec3& M2_R2,
						 Vec3& F2_I,
						 Vec3& M2_I,
						 Vec3& F1_I,
						 Vec3& M1_I,
						 doublereal& eps,
						 doublereal& eps_dot,
						 doublereal& delta,
						 doublereal& SoD,
						 doublereal& SoV,
						 doublereal& my,
						 doublereal& beta,
                         doublereal r,
                         doublereal alpha) const;
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	unsigned int iGetNumPrivData(void) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		      const VectorHandler& XCurr);
   	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
  private:
	const StructNode* m_pShaft;
	const StructNode* m_pBearing;
	Vec3 m_o1_R1;
	Vec3 m_o2_R2;
    bearing_data m_bdat;
	DriveOwner m_InitialAssemblyFactor;
    integer m_iNumGaussPoints;
	const doublereal* m_r;
	const doublereal* m_alpha;
    struct OutputOpt {
        doublereal r;
        doublereal alpha;
    } m_output[6];
    integer m_iNumOutputPoints;
	bool m_lambda;
	static const doublereal s_r1[1], s_alpha1[1];
	static const doublereal s_r2[2], s_alpha2[2];
	static const doublereal s_r3[3], s_alpha3[3];
	static const doublereal s_r6[6], s_alpha6[6];
};

// GAUSS-LEGENDRE integration points and weight factors according to K.J. Bathe 2002 table 5.6
// times 0.5 in order to account for an integration over the interval -0.5*b ... 0.5*b

const doublereal HydrodynamicPlainBearing::s_r1[1]     = {0.5 * 0.};
const doublereal HydrodynamicPlainBearing::s_alpha1[1] = {0.5 * 2};

const doublereal HydrodynamicPlainBearing::s_r2[2] = {0.5 * -0.577350269189626,
																	0.5 *  0.577350269189626};
const doublereal HydrodynamicPlainBearing::s_alpha2[2] = {0.5 * 1.,
																		0.5 * 1.};

const doublereal HydrodynamicPlainBearing::s_r3[3] = {0.5 * -0.774596669241483,
																	0.5 *  0.,
																	0.5 *  0.774596669241483};
const doublereal HydrodynamicPlainBearing::s_alpha3[3] = {0.5 * 0.555555555555556,
																		0.5 * 0.888888888888889,
																		0.5 * 0.555555555555556};

const doublereal HydrodynamicPlainBearing::s_r6[6] = {0.5 * -0.932469514203152,
																	0.5 * -0.661209386466265,
																	0.5 * -0.238619186083197,
																	0.5 *  0.238619186083197,
																	0.5 *  0.661209386466265,
																	0.5 *  0.932469514203152};
const doublereal HydrodynamicPlainBearing::s_alpha6[6] = {0.5 * 0.171324492379170,
																		0.5 * 0.360761573048139,
																		0.5 * 0.467913934572691,
																		0.5 * 0.467913934572691,
																		0.5 * 0.360761573048139,
																		0.5 * 0.171324492379170};

HydrodynamicPlainBearing::HydrodynamicPlainBearing(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: 	Elem(uLabel, flag(0)),
	UserDefinedElem(uLabel, pDO),
	m_pShaft(0),
	m_pBearing(0),
	m_o1_R1(0.,0.,0.),
	m_o2_R2(0.,0.,0.),
	m_InitialAssemblyFactor(0),
	m_iNumGaussPoints(0),
	m_r(0),
	m_alpha(0),
        m_iNumOutputPoints(0),
	m_lambda(true)
{
    std::memset(&m_bdat, 0, sizeof(m_bdat));
    std::memset(&m_output, 0, sizeof(m_output));
    
	// help
	if (HP.IsKeyWord("help"))
	{
		silent_cout(
			"\n"
			"Module: 	hydrodynamic_plain_bearing_with_offset\n"
			"\n"
			"	This module implements a hydrodynamic plain bearing according to\n"
			"\n"
			"   Hans Juergen Butenschoen 1976\n"
			"   Das hydrodynamische zylindrische Gleitlager\n"
			"	endlicher Breite unter instationaerer Belastung\n"
			"\n"
			"	hydrodynamic_plain_bearing_with_offset,\n"
			"		shaft, (label) <shaft node>,\n"
			"		[offset, (Vec3) <o1>,]\n"
			"		bearing, (label) <bearing node>,\n"
			"		[offset, (Vec3) <o2>,]\n"
			"		bearing_width, (real) <b>,\n"
			"		{shaft diameter, (real) <d> | bearing_diameter, (real) <D>,}\n"
			"		relative_clearance, (real) <Psi>,\n"
			"		oil_viscosity, (real) <eta>,\n"
			"		initial_assembly_factor, (DriveCaller),\n"
			"       [number of gauss points, <num_gauss_points>]\n"
            "       [output points, (integer) <num_output_points> [, {gauss point, (integer) <index_1> | custom, r, (real) <r_1>, alpha, (real) <alpha_1>}]]\n"
            "       [extend shaft to bearing center, {yes | no | (bool) <extend>}]\n"
			"\n"
			"   b ... bearing width [m]\n"
			"	d ... bearing diameter [m]\n"
			"	Psi ... relative clearance Psi = ( D - d ) / D [m/m]\n"
			"	eta ... dynamic oil viscosity [Pa*s]\n"
			"\n"
			<< std::endl);

		if (!HP.IsArg())
		{
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!HP.IsKeyWord("shaft")) 
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"shaft\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_pShaft = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	if (!m_pShaft) {
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if ( HP.IsKeyWord("offset") )
		m_o1_R1 = HP.GetPosRel(ReferenceFrame(m_pShaft));

	if ( !HP.IsKeyWord("bearing"))
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"bearing\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_pBearing = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	if (!m_pBearing) {
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if ( HP.IsKeyWord("offset") )
		m_o2_R2 = HP.GetPosRel(ReferenceFrame(m_pBearing));

	if ( !( HP.IsKeyWord("bearing_width") || HP.IsKeyWord("bearing" "width") ) )
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"bearing width\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}		
        
    m_bdat.b = HP.GetReal();

	bool bHaveD = false;

	if (HP.IsKeyWord("shaft" "diameter"))
	{
        m_bdat.d = HP.GetReal();
	}
	else
	{
		if ( !( HP.IsKeyWord("bearing_diameter") || HP.IsKeyWord("bearing" "diameter") ) )
		{
			silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"bearing diameter\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

        m_bdat.d = HP.GetReal();
		bHaveD = true;
	}

	if ( !( HP.IsKeyWord("relative_clearance") || HP.IsKeyWord("relative" "clearance") ) )
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"relative clearance\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
        
    m_bdat.Psi = HP.GetReal();

	if (bHaveD)
	{
        m_bdat.d *= (1 - m_bdat.Psi);
	}

	if ( !( HP.IsKeyWord("oil_viscosity") || HP.IsKeyWord("oil" "viscosity") ) )
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"oil viscosity\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
        
    m_bdat.eta = HP.GetReal();

	m_InitialAssemblyFactor.Set(( HP.IsKeyWord("initial_assembly_factor") || HP.IsKeyWord("initial" "assembly" "factor") ) ? HP.GetDriveCaller() : new OneDriveCaller());

    m_iNumGaussPoints = HP.IsKeyWord("number" "of" "gauss" "points") ? HP.GetInt() : 1;

#define CASE_GAUSS_POINT_NUM_(num)                                      \
    case num:                                                           \
        assert(m_iNumGaussPoints == sizeof(s_r##num)/sizeof(s_r##num[0])); \
        assert(m_iNumGaussPoints == sizeof(s_alpha##num)/sizeof(s_alpha##num[0])); \
        m_r = s_r##num;                                                 \
        m_alpha = s_alpha##num;                                         \
	break

    switch (m_iNumGaussPoints)
		{
	CASE_GAUSS_POINT_NUM_(1);
	CASE_GAUSS_POINT_NUM_(2);
	CASE_GAUSS_POINT_NUM_(3);
	CASE_GAUSS_POINT_NUM_(6);
    default:
        silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel()
                    << "): integration rule with " << m_iNumGaussPoints
                    << " gauss points not supported at line "
                    << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

    hydrodynamic_plain_bearing_init(m_bdat);

#undef CASE_GAUSS_POINT_NUM_
        
    ASSERT(m_iNumGaussPoints <= sizeof(m_output) / sizeof(m_output[0]));
        
    if (HP.IsKeyWord("output" "points")) {
        m_iNumOutputPoints = HP.GetInt();

        if (m_iNumOutputPoints < 0 || m_iNumOutputPoints > m_iNumGaussPoints) {
            silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel()
                        << "): number of points for output is not in range [0:" << m_iNumGaussPoints
                        << "] at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

        for (integer i = 0; i < m_iNumOutputPoints; ++i) {
            if (HP.IsKeyWord("gauss" "point")) {
                integer iGaussPoint = HP.GetInt();

                if (iGaussPoint < 1 || iGaussPoint > m_iNumGaussPoints) {
                    silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel()
                                << "): Gauss point index " << iGaussPoint
                                << " is not in range [1:" << m_iNumGaussPoints << "] at line "
                                << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

                m_output[i].r = m_r[iGaussPoint - 1];
                m_output[i].alpha = m_alpha[iGaussPoint - 1];
            } else if (HP.IsKeyWord("custom")) {
                if (!HP.IsKeyWord("r")) {
                    silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel()
                                << "): keyword \"r\" expected at line "
                                << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

                m_output[i].r = HP.GetReal();

                if (fabs(m_output[i].r) > 0.5) {
                    silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel()
                                << "): abs(r) = " << fabs(m_output[i].r)
                                << " > 0.5 is outside the bearing width at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (!HP.IsKeyWord("alpha")) {
                    silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel()
                                << "): keyword \"alpha\" expected at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

                m_output[i].alpha = HP.GetReal();

                if (m_output[i].alpha < 0 || m_output[i].alpha > 1) {
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel()
                                << "): alpha = " << m_output[i].alpha << " is not in range [0:1] at line "
					<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
            }
        }
    } else {
        m_iNumOutputPoints = m_iNumGaussPoints;

        for (integer i = 0; i < m_iNumOutputPoints; ++i) {
            m_output[i].r = m_r[i];
            m_output[i].alpha = m_alpha[i];
        }
    }

    if (HP.IsKeyWord("extend" "shaft" "to" "bearing" "center")) {
        m_lambda = HP.GetYesNoOrBool();
    }

    if (HP.IsKeyWord("epsilon" "max")) {
        m_bdat.eps_max = HP.GetReal();
    } else {
        m_bdat.eps_max = 0.999; // According to Butenschoen this model is valid until epsilon = 0.999
    }

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	std::ostream& out = pDM->GetLogFile();

	out << "hydrodynamic_plain_bearing_with_offset: "
		<< uLabel << " "
		<< m_pShaft->GetLabel() << " "
		<< m_o1_R1 << " "
		<< m_pBearing->GetLabel() << " "
		<< m_o2_R2 << " "
        << m_bdat.b << " "
        << m_bdat.d << " "
        << m_bdat.Psi << " "
        << m_bdat.eta << " ";

    out << m_iNumGaussPoints << " ";

    for (integer i = 0; i < m_iNumGaussPoints; ++i)
	{
        out << m_r[i] << " " << m_alpha[i] << " ";
	}

    out << m_iNumOutputPoints << " ";

    for (integer i = 0; i < m_iNumOutputPoints; ++i)
	{
        out << m_output[i].r << " " << m_output[i].alpha << " ";
	}

	out << std::endl;
}

HydrodynamicPlainBearing::~HydrodynamicPlainBearing(void)
{
    // destroy private dataman
}

void
HydrodynamicPlainBearing::Output(OutputHandler& OH) const
{
	if ( bToBeOutput() && OH.UseText(OutputHandler::LOADABLE) )
	{
		std::ostream& os = OH.Loadable();

		os << std::setw(8) << GetLabel() << ' '; 		// 0

        for (integer i = 0; i < m_iNumOutputPoints; ++i)
		{
			doublereal eps, eps_dot, delta, SoD, SoV, my, beta;

			Vec3 e_R2, e_dot_R2;
			doublereal omega_proj[2];
			Vec3 F2_R2, M2_R2;
			Vec3 F2_I, M2_I, F1_I, M1_I;

            ComputeResidual(e_R2, e_dot_R2, omega_proj, F2_R2, M2_R2, F2_I, M2_I, F1_I, M1_I, eps, eps_dot, delta, SoD, SoV, my, beta, m_output[i].r, m_output[i].alpha);

			// output the current state: the column layout is as follows
			// (all quantities are refered to the reference frame of the bearing node)

			// 1     2     3         4         5         6         7    8        9      10  11  12  13   14   15  16
			// e(1)  e(2)  e_dot(1)  e_dot(2)  omega(1)  omega(2)  eps  eps_dot  delta  Fx  Fy  Mz  SoD  SoV  my  beta

			//  0	label of the element
			//  1	absolute eccentricity in of the shaft in x direction
			//  2	absolute eccentricity in of the shaft in y direction
			//  3	difference of the absolute velocity between the shaft and the bearing in x direction
			//  4	difference of the absolute velocity between the shaft and the bearing in y direction
			//  5	absolute angular velocity of the shaft around the z axis
			//  6	absolute angular velocity of the bearing around the z axis
			//  7	relative eccentricity of the shaft
			//  8	time derivative of the relative eccentricity
			//  9	angle of minimum clearance
			// 10	force at the bearing in x direction
			// 11	force at the bearing in y direction
			// 12	torque at the bearing around the z axis
			// 13	Sommerfeld number for rotation
			// 14	Sommerfeld number for displacement	
			// 15	friction coefficient
			// 16	angle between minimum clearance and force of rotation

			os 	<< e_R2(1) << ' ' 			// 1
				<< e_R2(2) << ' ' 			// 2
				<< e_dot_R2(1) << ' '		// 3
				<< e_dot_R2(2) << ' '		// 4
				<< omega_proj[0] << ' '		// 5
				<< omega_proj[1] << ' '		// 6
				<< eps << ' ' 				// 7
				<< eps_dot << ' '			// 8
				<< delta << ' ' 			// 9
				<< F2_R2(1) << ' ' 			// 10
				<< F2_R2(2) << ' ' 			// 11
				<< M2_R2(3) << ' ' 			// 12
				<< SoD << ' ' 				// 13
				<< SoV << ' ' 				// 14
				<< my << ' ' 				// 15
				<< beta << ' ';				// 16
		}

		os << std::endl;
	}
}

void
HydrodynamicPlainBearing::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 12;
}

VariableSubMatrixHandler&
HydrodynamicPlainBearing::AssJac(VariableSubMatrixHandler& WorkMatV,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	integer iNumRows = 0;
	integer iNumCols = 0;

	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WorkMat = WorkMatV.SetFull();

	WorkMat.ResizeReset(iNumRows, iNumCols);

	integer intShaftPositionIndex = m_pShaft->iGetFirstPositionIndex();
	integer intShaftMomentumIndex = m_pShaft->iGetFirstMomentumIndex();
	integer intBearingPositionIndex = m_pBearing->iGetFirstPositionIndex();
	integer intBearingMomentumIndex = m_pBearing->iGetFirstMomentumIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++)
	{
		WorkMat.PutRowIndex(iCnt, intShaftMomentumIndex + iCnt);
		WorkMat.PutColIndex(iCnt, intShaftPositionIndex + iCnt);
		WorkMat.PutRowIndex(iCnt+6,intBearingMomentumIndex + iCnt);
		WorkMat.PutColIndex(iCnt+6,intBearingPositionIndex + iCnt);
	}

	const Vec3& X1 = m_pShaft->GetXCurr();
	const Vec3& X2 = m_pBearing->GetXCurr();
	const Vec3& X1_dot = m_pShaft->GetVCurr();
	const Vec3& X2_dot = m_pBearing->GetVCurr();
	const Mat3x3& R1 = m_pShaft->GetRCurr();
	const Mat3x3& R1_0 = m_pShaft->GetRRef();
	const Mat3x3& R2 = m_pBearing->GetRCurr();
	const Mat3x3& R2_0 = m_pBearing->GetRRef();
	const Vec3& omega1 = m_pShaft->GetWCurr();
	const Vec3& omega1_ref = m_pShaft->GetWRef();
	const Vec3& omega2 = m_pBearing->GetWCurr();
	const Vec3& omega2_ref = m_pBearing->GetWRef();

    for (integer i = 0; i < m_iNumGaussPoints; ++i)
	{
		Vec3 o1_R1 = m_o1_R1;
        o1_R1(3) += m_r[i] * m_bdat.b;
        Vec3 o2_R2 = m_o2_R2;
        o2_R2(3) += m_r[i] * m_bdat.b;
        const Vec3 v_R2 = R2.MulTV( X1 - X2 + R1 * o1_R1 ) - o2_R2;
		const Vec3 d1_R2 = R2.MulTV(R1.GetCol(3));

		const doublereal lambda = m_lambda ? -v_R2(3) / d1_R2(3) : 0.;

		Vec3 e_R2 = v_R2 + d1_R2 * lambda;

		// FIXME: nan values if e_R2(1) == 0 && e_R2(2) == 0
		if ( e_R2(1) == 0.0 && e_R2(2) == 0.0 )
            e_R2(1) = std::numeric_limits<doublereal>::epsilon() * m_bdat.d * m_bdat.Psi / 2.;

		const Vec3 v_dot_R2 = R2.MulTV( X1_dot - X2_dot + omega2.Cross(X2 - X1) + ( omega1 - omega2 ).Cross( R1 * o1_R1 ) );
		const Vec3 d1_dot_R2 = R2.MulTV( ( omega1 - omega2 ).Cross( R1.GetCol(3) ) );
		const doublereal lambda_dot = m_lambda ? -v_dot_R2(3) / d1_R2(3) + v_R2(3) / std::pow(d1_R2(3),2) * d1_dot_R2(3) : 0.;

        Vec3 e_dot_R2 = R2.MulTV( X1_dot - X2_dot + omega1.Cross( R1 * o1_R1 ) - omega2.Cross( R2 * o2_R2 )
								+ R1.GetCol(3) * lambda_dot + omega1.Cross(R1.GetCol(3)) * lambda);

		// FIXME: nan values if e_dot_R2(1) == 0 && e_dot_R2(2) == 0
		if ( e_dot_R2(1) == 0.0 && e_dot_R2(2) == 0.0 )
            e_dot_R2(1) = std::numeric_limits<doublereal>::epsilon() * m_bdat.d * m_bdat.Psi / 2.;

		const Vec3 lambda_d1_R1 = Vec3(0.,0.,lambda);
		const Vec3 l1_R1 = o1_R1 + lambda_d1_R1;
		const Vec3 l1_I = R1 * l1_R1;
        const Vec3 l2_R2 = o2_R2 + e_R2;
		const doublereal omega_proj[2] = { R2.GetCol(3).Dot(omega1),
						   R2.GetCol(3).Dot(omega2) };

		const doublereal e_[2] = { e_R2(1), e_R2(2) }, e_dot_[2] = { e_dot_R2(1), e_dot_R2(2) };

		//												  |	0   1   2   3   4   5  |
		static const doublereal ed[2][NBDIRSMAX] 	 	  = { { 1., 0., 0., 0., 0., 0. },	// 0 | inner derivative of the eccentricity of the shaft in direction 1 of the reference frame of the bearing (R2)
														  { 0., 1., 0., 0., 0., 0. } };	// 1 | inner derivative of the eccentricity of the shaft in direction 2 of the reference frame of the bearing (R2)
		static const doublereal e_dotd[2][NBDIRSMAX] 	  = { { 0., 0., 1., 0., 0., 0. },	// 0 | inner derivative of relative velocity of the shaft in direction 1 of the reference frame of the bearing (R2)
														  { 0., 0., 0., 1., 0., 0. } }; // 1 | inner derivative of relative velocity of the shaft in direction 2 of the reference frame of the bearing (R2)
		static const doublereal omega_projd[2][NBDIRSMAX] = { { 0., 0., 0., 0., 1., 0. },	// 0 | inner derivative of the angular velocity of the shaft in direction 3 of the reference frame of the bearing (R2)
														  { 0., 0., 0., 0., 0., 1. } }; // 1 | inner derivative of the angular velocity of the bearing in direction 3 of the reference frame of the bearing (R2)

		doublereal k[3];			// force vector at the bearing (not used for the evaluation of the jacobian matrix)
		doublereal kd[3][NBDIRSMAX];	// variation of the force vector at the bearing with respect to the eccentricity of the shaft e and the relative velocity of the shaft

		// kd = { { diff(k[0],e[0]), diff(k[0],e[1]), diff(k[0],e_dot[0]), diff(k[0],e_dot[1]), diff(k[0],omega_proj[0]), diff(k[0],omega_proj[1]) },
		//        { diff(k[1],e[0]), diff(k[1],e[1]), diff(k[1],e_dot[0]), diff(k[1],e_dot[1]), diff(k[1],omega_proj[0]), diff(k[1],omega_proj[1]) },
		//        { diff(k[2],e[0]), diff(k[2],e[1]), diff(k[2],e_dot[0]), diff(k[2],e_dot[1]), diff(k[2],omega_proj[0]), diff(k[2],omega_proj[1]) } };

		doublereal eps, eps_dot, delta, SoD, SoV, my, beta;

        hydrodynamic_plain_bearing_force_dv(m_bdat, omega_proj, omega_projd, e_, ed, e_dot_, e_dotd, k, kd, eps, eps_dot, delta, SoD, SoV, my, beta, NBDIRSMAX);

		Vec3 F2_R2;					// F2^(R2) = f(e^(R2), diff(e^(R2),t), omega1_proj, omega2_proj)

		F2_R2(1) = k[0];
		F2_R2(2) = k[1];
		F2_R2(3) = 0.;

		Mat3x3 dF2_R2_de_R2;		// diff(F2^(R2), e^(R2))

		dF2_R2_de_R2(1,1) = kd[0][0]; dF2_R2_de_R2(1,2) = kd[0][1]; dF2_R2_de_R2(1,3) = 0.;
		dF2_R2_de_R2(2,1) = kd[1][0]; dF2_R2_de_R2(2,2) = kd[1][1]; dF2_R2_de_R2(2,3) = 0.;
		dF2_R2_de_R2(3,1) = 	  0.; dF2_R2_de_R2(3,2) = 		0.; dF2_R2_de_R2(3,3) = 0.;

		Mat3x3 dF2_R2_de_dot_R2;	// diff(F2^(R2), diff(e^(R2),t))

		dF2_R2_de_dot_R2(1,1) = kd[0][2]; dF2_R2_de_dot_R2(1,2) = kd[0][3]; dF2_R2_de_dot_R2(1,3) = 0.;
		dF2_R2_de_dot_R2(2,1) = kd[1][2]; dF2_R2_de_dot_R2(2,2) = kd[1][3]; dF2_R2_de_dot_R2(2,3) = 0.;
		dF2_R2_de_dot_R2(3,1) =		  0.; dF2_R2_de_dot_R2(3,2) =		 0; dF2_R2_de_dot_R2(3,3) = 0.;

		Vec3 dF2_R2_domega1_proj;	// diff(F2^(R2), omega1_proj)

		dF2_R2_domega1_proj(1) = kd[0][4];
		dF2_R2_domega1_proj(2) = kd[1][4];
		dF2_R2_domega1_proj(3) = 	   0.;

		Vec3 dF2_R2_domega2_proj;	// diff(F2^(R2), omega2_proj)

		dF2_R2_domega2_proj(1) = kd[0][5];
		dF2_R2_domega2_proj(2) = kd[1][5];
		dF2_R2_domega2_proj(3) =       0.;

		Vec3 M2_R2;					// M2^(R2) = f(e^(R2), diff(e^(R2),t), omega1_proj, omega2_proj)

		M2_R2(1) = 0.;
		M2_R2(2) = 0.;
		M2_R2(3) = k[2];

		Mat3x3 dM2_R2_de_R2;		// diff(M2^(R2), e^(R2)

		dM2_R2_de_R2(1,1) = 	  0.; 	dM2_R2_de_R2(1,2) = 	  0.; 	dM2_R2_de_R2(1,3) = 0.;
		dM2_R2_de_R2(2,1) = 	  0.; 	dM2_R2_de_R2(2,2) = 	  0.; 	dM2_R2_de_R2(2,3) = 0.;
		dM2_R2_de_R2(3,1) = kd[2][0]; 	dM2_R2_de_R2(3,2) = kd[2][1]; 	dM2_R2_de_R2(3,3) = 0.;

		Mat3x3 dM2_R2_de_dot_R2;	// diff(M2^(R2), diff(e^(R2),t)

		dM2_R2_de_dot_R2(1,1) = 	  0.; 	dM2_R2_de_dot_R2(1,2) = 	  0.; 	dM2_R2_de_dot_R2(1,3) = 0.;
		dM2_R2_de_dot_R2(2,1) = 	  0.; 	dM2_R2_de_dot_R2(2,2) = 	  0.; 	dM2_R2_de_dot_R2(2,3) = 0.;
		dM2_R2_de_dot_R2(3,1) = kd[2][2]; 	dM2_R2_de_dot_R2(3,2) = kd[2][3]; 	dM2_R2_de_dot_R2(3,3) = 0.;

		Vec3 dM2_R2_domega1_proj;	// diff(M2^(R2), omega1_proj)

		dM2_R2_domega1_proj(1) = 0.;
		dM2_R2_domega1_proj(2) = 0.;
		dM2_R2_domega1_proj(3) = kd[2][4];

		Vec3 dM2_R2_domega2_proj;	// diff(M2^(R2), omega2_proj)

		dM2_R2_domega2_proj(1) = 0.;
		dM2_R2_domega2_proj(2) = 0.;
		dM2_R2_domega2_proj(3) = kd[2][5];

		const doublereal alpha = m_alpha[i] * m_InitialAssemblyFactor.dGet();

		F2_R2 *= alpha;
		dF2_R2_de_R2 *= alpha;
		dF2_R2_de_dot_R2 *= alpha;
		dF2_R2_domega1_proj *= alpha;
		dF2_R2_domega2_proj *= alpha;

		M2_R2 *= alpha;
		dM2_R2_de_R2 *= alpha;
		dM2_R2_de_dot_R2 *= alpha;
		dM2_R2_domega1_proj *= alpha;
		dM2_R2_domega2_proj *= alpha;

		const Mat3x3 R2_T = R2.Transpose();

		const Vec3 F2_I = R2 * F2_R2;
		const Vec3 M2_I = R2 * ( l2_R2.Cross( F2_R2 ) + M2_R2 );
		const Vec3 F1_I = -F2_I;
		const Vec3 M1_I = -l1_I.Cross( F2_I ) - R2 * M2_R2;

		const Mat3x3 dv_dot_R2_dX1 = -R2.MulTM(Mat3x3(MatCross, omega2));	// diff(diff(v^(R2),t),X1)
		const Mat3x3& dv_R2_dX1 = R2_T;				// diff(v^(R2),X1)
		const Vec3 dlambda_dot_dX1 = m_lambda ? -dv_dot_R2_dX1.GetRow(3) / d1_R2(3) + dv_R2_dX1.GetRow(3) / std::pow(d1_R2(3),2) * d1_dot_R2(3) : Zero3; // diff(diff(lambda,t),X1)

		const Vec3 dlambda_dX1 = m_lambda ? -dv_R2_dX1.GetRow(3) / d1_R2(3) : Zero3;

		const Mat3x3 de_dot_R2_dX1 = R2.MulTM( R1.GetCol(3).Tens(dlambda_dot_dX1) + omega1.Cross( R1.GetCol(3) ).Tens(dlambda_dX1) );

		const Mat3x3 dv_dot_R2_domega1 = -R2.MulTM(Mat3x3(MatCross, R1 * o1_R1));
		const Mat3x3 domega1_dg1 = -Mat3x3(MatCross, omega1_ref);
		const Mat3x3 dv_R2_dg1 = -R2.MulTM( Mat3x3(MatCross, R1_0 * o1_R1) );
		const Mat3x3 dv_dot_R2_dg1 = dv_dot_R2_domega1 * domega1_dg1 - R2.MulTM( ( omega1 - omega2 ).Cross(Mat3x3(MatCross, R1_0 * o1_R1)) );
		const Mat3x3 dd1_R2_dg1 = -R2.MulTM(Mat3x3(MatCross, R1_0.GetCol(3)));
		const Mat3x3 dd1_dot_R2_dg1 = R2.MulTM( R1.GetCol(3).Cross(Mat3x3(MatCross, omega1_ref)) - ( omega1 - omega2 ).Cross(Mat3x3(MatCross, R1_0.GetCol(3))) );
		const Vec3 dlambda_dot_dg1 = m_lambda ? -dv_dot_R2_dg1.GetRow(3) / d1_R2(3)
						+ dd1_R2_dg1.GetRow(3) * (v_dot_R2(3) / std::pow(d1_R2(3), 2))
						+ dv_R2_dg1.GetRow(3) / std::pow(d1_R2(3), 2) * d1_dot_R2(3)
						- dd1_R2_dg1.GetRow(3) * (2. * v_R2(3) / std::pow(d1_R2(3), 3) * d1_dot_R2(3))
						+ dd1_dot_R2_dg1.GetRow(3) * (v_R2(3) / std::pow(d1_R2(3), 2)) : Zero3;
		const Vec3 dlambda_dg1 = m_lambda ? -dv_R2_dg1.GetRow(3) / d1_R2(3) + dd1_R2_dg1.GetRow(3) * (v_R2(3) / std::pow( d1_R2(3), 2 )) : Zero3;

		const Mat3x3 de_dot_R2_dg1 = R2.MulTM( -omega1.Cross(Mat3x3( MatCross, R1_0 * o1_R1 )) + ( R1 * o1_R1 ).Cross(Mat3x3(MatCross, omega1_ref))
									+ R1.GetCol(3).Tens(dlambda_dot_dg1) - Mat3x3( MatCross, R1_0.GetCol(3) ) * lambda_dot
									+ omega1.Cross( R1.GetCol(3) ).Tens(dlambda_dg1) + R1.GetCol(3).Cross(Mat3x3(MatCross, omega1_ref)) * lambda
									- omega1.Cross(Mat3x3( MatCross, R1_0.GetCol(3) )) * lambda );
		const Mat3x3 dv_dot_R2_dX2 = R2.MulTM(Mat3x3(MatCross, omega2));
		const Mat3x3 dv_R2_dX2 = -R2_T;
		const Vec3 dlambda_dot_dX2 = m_lambda ? -dv_dot_R2_dX2.GetRow(3) / d1_R2(3) + dv_R2_dX2.GetRow(3) / std::pow(d1_R2(3),2) * d1_dot_R2(3) : Zero3;

		const Vec3 dlambda_dX2 = m_lambda ? -dv_R2_dX2.GetRow(3) / d1_R2(3) : Zero3;
		const Mat3x3 de_dot_R2_dX2 = R2.MulTM( R1.GetCol(3).Tens(dlambda_dot_dX2) + omega1.Cross( R1.GetCol(3) ).Tens(dlambda_dX2) );
		const Mat3x3 dv_R2_dg2 = R2_0.MulTM( Mat3x3( MatCross, X1 - X2 + R1 * o1_R1 ) );
		const Mat3x3 dd1_R2_dg2 = R2_0.MulTM( Mat3x3( MatCross, R1 * o1_R1 ) );
		const Mat3x3 dv_dot_R2_dg2 = R2_0.MulTM( Mat3x3( MatCross, X1_dot - X2_dot + omega2.Cross( X2 - X1 ) + ( omega1 - omega2 ).Cross( R1 * o1_R1 ) ) )
					   + R2.MulTM( ( X2 - X1 - R1 * o1_R1 ).Cross(Mat3x3(MatCross, omega2_ref)) );
		const Vec3 dlambda_dg2 = m_lambda ? -dv_R2_dg2.GetRow(3) / d1_R2(3) + dd1_R2_dg2.GetRow(3) * (v_R2(3) / std::pow(d1_R2(3),2)) : Zero3;

		const Mat3x3 dd1_dot_R2_dg2 = R2_0.MulTM( Mat3x3( MatCross, ( omega1 - omega2 ).Cross( R1.GetCol(3) ) ) )
									- R2.MulTM( Mat3x3( MatCross, R1.GetCol(3) ) * Mat3x3(MatCross, omega2_ref) );
		const Vec3 dlambda_dot_dg2 = m_lambda ? -dv_dot_R2_dg2.GetRow(3) / d1_R2(3) + dd1_R2_dg2.GetRow(3) * (v_dot_R2(3) / std::pow(d1_R2(3), 2))
						+ dv_R2_dg2.GetRow(3) / std::pow(d1_R2(3),2) * d1_dot_R2(3)
						- dd1_R2_dg2.GetRow(3) * (2. * v_R2(3) / std::pow( d1_R2(3), 3) * d1_dot_R2(3))
						+ dd1_dot_R2_dg2.GetRow(3) * (v_R2(3) / std::pow(d1_R2(3),2)) : Zero3;

        const Mat3x3 de_dot_R2_dg2 = R2_0.MulTM( Mat3x3( MatCross, X1_dot - X2_dot + omega1.Cross(R1 * o1_R1) - omega2.Cross( R2 * o2_R2 )
									+ R1.GetCol(3) * lambda_dot + omega1.Cross(R1.GetCol(3)) * lambda))
            + R2.MulTM( -( R2 * o2_R2 ).Cross(Mat3x3(MatCross, omega2_ref)) + omega2.Cross(Mat3x3( MatCross, R2_0 * o2_R2 ))
									+ R1.GetCol(3).Tens(dlambda_dot_dg2) + omega1.Cross(R1.GetCol(3)).Tens(dlambda_dg2));
		const Mat3x3& dv_dot_R2_dX1_dot = R2_T;
		const Vec3 dlambda_dot_dX1_dot = m_lambda ? -dv_dot_R2_dX1_dot.GetRow(3) / d1_R2(3) : Zero3;

		const Mat3x3 de_dot_R2_dX1_dot = R2.MulTM( Eye3 + R1.GetCol(3).Tens(dlambda_dot_dX1_dot));

		const Mat3x3 dv_dot_R2_dg1_dot = -R2.MulTM( Mat3x3( MatCross, R1 * o1_R1 ) );
		const Mat3x3 dd1_dot_R2_dg1_dot = -R2.MulTM( Mat3x3( MatCross, R1.GetCol(3) ) );
		const Vec3 dlambda_dot_dg1_dot = m_lambda ? -dv_dot_R2_dg1_dot.GetRow(3) / d1_R2(3) + dd1_dot_R2_dg1_dot.GetRow(3) * (v_R2(3) / std::pow(d1_R2(3),2)) : Zero3;

		const Mat3x3 de_dot_R2_dg1_dot = R2.MulTM( -Mat3x3( MatCross, R1 * o1_R1 ) + R1.GetCol(3).Tens(dlambda_dot_dg1_dot) - Mat3x3(MatCross, R1.GetCol(3) * lambda) );

		const Mat3x3 dv_dot_dX2_dot = -R2_T;
		const Vec3 dlambda_dot_dX2_dot = m_lambda ? -dv_dot_dX2_dot.GetRow(3) / d1_R2(3) : Zero3;

		const Mat3x3 de_dot_R2_dX2_dot = R2.MulTM( -Eye3 + R1.GetCol(3).Tens(dlambda_dot_dX2_dot));

		const Mat3x3 dv_dot_R2_dg2_dot = R2.MulTM( Mat3x3( MatCross, R1 * o1_R1 + X1 - X2 ));
		const Mat3x3 dd1_dot_R2_dg2_dot = R2.MulTM(Mat3x3(MatCross, R1.GetCol(3)));
		const Vec3 dlambda_dot_dg2_dot = m_lambda ? -dv_dot_R2_dg2_dot.GetRow(3) / d1_R2(3) + dd1_dot_R2_dg2_dot.GetRow(3) * (v_R2(3) / std::pow(d1_R2(3),2)) : Zero3;

        const Mat3x3 de_dot_R2_dg2_dot = R2.MulTM( Mat3x3( MatCross, R2 * o2_R2 ) + R1.GetCol(3).Tens( dlambda_dot_dg2_dot) );

		const Mat3x3 de_R2_dX1 = dv_R2_dX1 + d1_R2.Tens( dlambda_dX1 );

		const Mat3x3 de_R2_dg1 = dv_R2_dg1 + d1_R2.Tens(dlambda_dg1) + dd1_R2_dg1 * lambda;

		const Mat3x3 de_R2_dX2 = dv_R2_dX2 + d1_R2.Tens(dlambda_dX2);

		const Mat3x3 de_R2_dg2 = dv_R2_dg2 + d1_R2.Tens( dlambda_dg2 ) + dd1_R2_dg2 * lambda;

		const Vec3 domega1_proj_dg1 = omega1_ref.Cross( R2.GetCol(3) );
		const Vec3 domega1_proj_dg1_dot = R2.GetCol(3);
		const Vec3 domega1_proj_dg2 = -omega1.Cross( R2_0.GetCol(3) );
		const Vec3 domega2_proj_dg2 = -omega2.Cross(R2_0.GetCol(3)) + omega2_ref.Cross(R2.GetCol(3));
		const Vec3 domega2_proj_dg2_dot = R2.GetCol(3);

		const Mat3x3 dF2_R2_dX1 = dF2_R2_de_R2 * de_R2_dX1 + dF2_R2_de_dot_R2 * de_dot_R2_dX1;

		const Mat3x3 dF2_R2_dg1 = dF2_R2_de_R2 * de_R2_dg1 + dF2_R2_de_dot_R2 * de_dot_R2_dg1
								+ dF2_R2_domega1_proj.Tens(domega1_proj_dg1);
		const Mat3x3 dF2_R2_dX2 = dF2_R2_de_R2 * de_R2_dX2 + dF2_R2_de_dot_R2 * de_dot_R2_dX2;
		const Mat3x3 dF2_R2_dg2 = dF2_R2_de_R2 * de_R2_dg2 + dF2_R2_de_dot_R2 * de_dot_R2_dg2
								+ dF2_R2_domega1_proj.Tens(domega1_proj_dg2) + dF2_R2_domega2_proj.Tens(domega2_proj_dg2);

		const Mat3x3 dF2_I_dX1 = R2 * dF2_R2_dX1;
		const Mat3x3 dF2_I_dg1 = R2 * dF2_R2_dg1;
		const Mat3x3 dF2_I_dX2 = R2 * dF2_R2_dX2;
		const Mat3x3 dF2_I_dg2 = -Mat3x3( MatCross, R2_0 * F2_R2 ) + R2 * dF2_R2_dg2;

		const Mat3x3 dF1_I_dX1 = -dF2_I_dX1;
		const Mat3x3 dF1_I_dg1 = -dF2_I_dg1;
		const Mat3x3 dF1_I_dX2 = -dF2_I_dX2;
		const Mat3x3 dF1_I_dg2 = -dF2_I_dg2;

		const Mat3x3 dM2_R2_dX1 = dM2_R2_de_R2 * de_R2_dX1 + dM2_R2_de_dot_R2 * de_dot_R2_dX1;
		const Mat3x3 dM2_I_dX1 = R2 * ( -F2_R2.Cross(de_R2_dX1) + l2_R2.Cross(dF2_R2_dX1) + dM2_R2_dX1 );

		const Mat3x3 dM2_R2_dg1 = dM2_R2_de_R2 * de_R2_dg1 + dM2_R2_de_dot_R2 * de_dot_R2_dg1 + dM2_R2_domega1_proj.Tens(domega1_proj_dg1);
		const Mat3x3 dM2_I_dg1 = R2 * ( -F2_R2.Cross(de_R2_dg1) + l2_R2.Cross(dF2_R2_dg1) + dM2_R2_dg1 );

		const Mat3x3 dM2_R2_dX2 = dM2_R2_de_R2 * de_R2_dX2 + dM2_R2_de_dot_R2 * de_dot_R2_dX2;
		const Mat3x3 dM2_I_dX2 = R2 * ( -F2_R2.Cross( de_R2_dX2 ) + l2_R2.Cross(dF2_R2_dX2) + dM2_R2_dX2 );

		const Mat3x3 dM2_R2_dg2 = dM2_R2_de_R2 * de_R2_dg2 + dM2_R2_de_dot_R2 * de_dot_R2_dg2
								+ dM2_R2_domega1_proj.Tens(domega1_proj_dg2) + dM2_R2_domega2_proj.Tens(domega2_proj_dg2);
		const Mat3x3 dM2_I_dg2 = -Mat3x3( MatCross, R2_0 * ( l2_R2.Cross(F2_R2) + M2_R2 ) )
								+ R2 * ( -F2_R2.Cross(de_R2_dg2) + l2_R2.Cross(dF2_R2_dg2) + dM2_R2_dg2 );

		const Mat3x3 dM1_I_dX1 = ( R2 * F2_R2 ).Cross( R1.GetCol(3) ).Tens( dlambda_dX1 )
							   - l1_I.Cross( R2 * dF2_R2_dX1 ) - R2 * dM2_R2_dX1;

		const Mat3x3 dM1_I_dg1 = ( R2 * F2_R2 ).Cross( R1.GetCol(3).Tens(dlambda_dg1) - Mat3x3( MatCross, R1_0 * l1_R1 ) )
								- l1_I.Cross( R2 * dF2_R2_dg1 ) - R2 * dM2_R2_dg1;

		const Mat3x3 dM1_I_dX2 = ( R2 * F2_R2 ).Cross( R1.GetCol(3).Tens(dlambda_dX2) )
							   - l1_I.Cross( R2 * dF2_R2_dX2 ) - R2 * dM2_R2_dX2;

		const Mat3x3 dM1_I_dg2 = ( R2 * F2_R2 ).Cross( R1.GetCol(3).Tens(dlambda_dg2) )
							   + l1_I.Cross( Mat3x3( MatCross, R2_0 * F2_R2 ) - R2 * dF2_R2_dg2 )
							   + Mat3x3( MatCross, R2_0 * M2_R2 ) - R2 * dM2_R2_dg2;

		const Mat3x3 dF2_R2_dX1_dot = dF2_R2_de_dot_R2 * de_dot_R2_dX1_dot;
		const Mat3x3 dF2_I_dX1_dot = R2 * dF2_R2_dX1_dot;

		const Mat3x3 dF2_R2_dg1_dot = dF2_R2_de_dot_R2 * de_dot_R2_dg1_dot + dF2_R2_domega1_proj.Tens( domega1_proj_dg1_dot );
		const Mat3x3 dF2_I_dg1_dot = R2 * dF2_R2_dg1_dot;

		const Mat3x3 dF2_R2_dX2_dot = dF2_R2_de_dot_R2 * de_dot_R2_dX2_dot;
		const Mat3x3 dF2_I_dX2_dot = R2 * dF2_R2_dX2_dot;

		const Mat3x3 dF2_R2_dg2_dot = dF2_R2_de_dot_R2 * de_dot_R2_dg2_dot + dF2_R2_domega2_proj.Tens( domega2_proj_dg2_dot );
		const Mat3x3 dF2_I_dg2_dot = R2 * dF2_R2_dg2_dot;

		const Mat3x3 dF1_I_dX1_dot = -dF2_I_dX1_dot;
		const Mat3x3 dF1_I_dg1_dot = -dF2_I_dg1_dot;
		const Mat3x3 dF1_I_dX2_dot = -dF2_I_dX2_dot;
		const Mat3x3 dF1_I_dg2_dot = -dF2_I_dg2_dot;

		const Mat3x3 dM2_R2_dX1_dot = dM2_R2_de_dot_R2 * de_dot_R2_dX1_dot;
		const Mat3x3 dM2_I_dX1_dot = R2 * ( l2_R2.Cross( dF2_R2_dX1_dot) + dM2_R2_dX1_dot );

		const Mat3x3 dM2_R2_dg1_dot = dM2_R2_de_dot_R2 * de_dot_R2_dg1_dot + dM2_R2_domega1_proj.Tens( domega1_proj_dg1_dot );
		const Mat3x3 dM2_I_dg1_dot = R2 * ( l2_R2.Cross( dF2_R2_dg1_dot ) + dM2_R2_dg1_dot );

		const Mat3x3 dM2_R2_dX2_dot = dM2_R2_de_dot_R2 * de_dot_R2_dX2_dot;
		const Mat3x3 dM2_I_dX2_dot = R2 * ( l2_R2.Cross(dF2_R2_dX2_dot) + dM2_R2_dX2_dot );

		const Mat3x3 dM2_R2_dg2_dot = dM2_R2_de_dot_R2 * de_dot_R2_dg2_dot + dM2_R2_domega2_proj.Tens(domega2_proj_dg2_dot);
		const Mat3x3 dM2_I_dg2_dot = R2 * ( l2_R2.Cross( dF2_R2_dg2_dot ) + dM2_R2_dg2_dot );

		const Mat3x3 dM1_I_dX1_dot = -l1_I.Cross( R2 * dF2_R2_dX1_dot ) - R2 * dM2_R2_dX1_dot;
		const Mat3x3 dM1_I_dg1_dot = -l1_I.Cross( R2 * dF2_R2_dg1_dot ) - R2 * dM2_R2_dg1_dot;
		const Mat3x3 dM1_I_dX2_dot = -l1_I.Cross( R2 * dF2_R2_dX2_dot ) - R2 * dM2_R2_dX2_dot;
		const Mat3x3 dM1_I_dg2_dot = -l1_I.Cross( R2 * dF2_R2_dg2_dot ) - R2 * dM2_R2_dg2_dot;

			/*
			 *                    1,           4,           7,          10                    1,       4,       7,      10
			 *        | dF1/dX1_dot, dF1/dg1_dot, dF1/dX2_dot, dF1/dg2_dot |          | dF1/dX1, dF1/dg1, dF1/dX2, dF1/dg2 |  1
			 *        | dM1/dX1_dot, dM1/dg1_dot, dM1/dX2_dot, dM1/dg2_dot |          | dM1/dX1, dM1/dg1, dM1/dX2, dM1/dg2 |  4
			 * Jac = -|                                                    | -dCoef * |                                    |
			 *        | dF2/dX1_dot, dF2/dg1_dot, dF2/dX2_dot, dF2/dg2_dot |          | dF2/dX1, dF2/dg1, dF2/dX2, dF2/dg2 |  7
			 *        | dM2/dX1_dot, dM2/dg1_dot, dM2/dX2_dot, dM2/dg2_dot |          | dM2/dX1, dM2/dg1, dM2/dX2, dM2/dg2 | 10
			*/
		WorkMat.Sub(  1,  1, dF1_I_dX1_dot + dF1_I_dX1 * dCoef );
		WorkMat.Sub(  1,  4, dF1_I_dg1_dot + dF1_I_dg1 * dCoef );
		WorkMat.Sub(  1,  7, dF1_I_dX2_dot + dF1_I_dX2 * dCoef );
		WorkMat.Sub(  1, 10, dF1_I_dg2_dot + dF1_I_dg2 * dCoef );

		WorkMat.Sub(  4,  1, dM1_I_dX1_dot + dM1_I_dX1 * dCoef );
		WorkMat.Sub(  4,  4, dM1_I_dg1_dot + dM1_I_dg1 * dCoef );
		WorkMat.Sub(  4,  7, dM1_I_dX2_dot + dM1_I_dX2 * dCoef );
		WorkMat.Sub(  4, 10, dM1_I_dg2_dot + dM1_I_dg2 * dCoef );

		WorkMat.Sub(  7,  1, dF2_I_dX1_dot + dF2_I_dX1 * dCoef );
		WorkMat.Sub(  7,  4, dF2_I_dg1_dot + dF2_I_dg1 * dCoef );
		WorkMat.Sub(  7,  7, dF2_I_dX2_dot + dF2_I_dX2 * dCoef );
		WorkMat.Sub(  7, 10, dF2_I_dg2_dot + dF2_I_dg2 * dCoef );

		WorkMat.Sub( 10,  1, dM2_I_dX1_dot + dM2_I_dX1 * dCoef );
		WorkMat.Sub( 10,  4, dM2_I_dg1_dot + dM2_I_dg1 * dCoef );
		WorkMat.Sub( 10,  7, dM2_I_dX2_dot + dM2_I_dX2 * dCoef );
		WorkMat.Sub( 10, 10, dM2_I_dg2_dot + dM2_I_dg2 * dCoef );
	}
#ifdef DEBUG
	std::cerr << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << ":" <<  "Jac=" << std::endl << WorkMat << std::endl;
#endif
	return WorkMatV;
}

SubVectorHandler&
HydrodynamicPlainBearing::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// resize residual
	integer iNumRows = 0;
	integer iNumCols = 0;

	WorkSpaceDim(&iNumRows, &iNumCols);

	WorkVec.ResizeReset(iNumRows);

	const integer intShaftMomentumIndex = m_pShaft->iGetFirstMomentumIndex();
	const integer intBearingMomentumIndex = m_pBearing->iGetFirstMomentumIndex();

	// equations indexes
	for ( int iCnt = 1; iCnt <= 6; ++iCnt)
	{
		WorkVec.PutRowIndex(iCnt,  intShaftMomentumIndex + iCnt);
		WorkVec.PutRowIndex(iCnt+6,intBearingMomentumIndex + iCnt);
	}

    for (integer i = 0; i < m_iNumGaussPoints; ++i)
	{
		doublereal eps, eps_dot, delta, SoD, SoV, my, beta;

		Vec3 e_R2, e_dot_R2;
		doublereal omega_proj[2];
		Vec3 F2_R2, M2_R2;
		Vec3 F2_I, M2_I, F1_I, M1_I;

        ComputeResidual(e_R2, e_dot_R2, omega_proj, F2_R2, M2_R2, F2_I, M2_I, F1_I, M1_I, eps, eps_dot, delta, SoD, SoV, my, beta, m_r[i], m_alpha[i]);
		// 1     2     3     4     5     6     7     8     9     10    11    12
		// F1(1) F1(2) F1(3) M1(1) M1(2) M1(3) F2(1) F2(2) F2(3) M2(1) M2(2) M2(3)
		WorkVec.Add(1,  F1_I);
		WorkVec.Add(4,  M1_I);
		WorkVec.Add(7,  F2_I);
		WorkVec.Add(10, M2_I);
	}

#ifdef DEBUG
	std::cerr << __FILE__ << ":" << __LINE__ << ":" << __PRETTY_FUNCTION__ << ": Res=" << std::endl;
	std::cerr << WorkVec << std::endl;
#endif
	return WorkVec;
}

void HydrodynamicPlainBearing::ComputeResidual(Vec3& e_R2, Vec3& e_dot_R2,doublereal omega_proj[2],Vec3& F2_R2,Vec3& M2_R2,Vec3& F2_I,Vec3& M2_I,Vec3& F1_I,Vec3& M1_I,doublereal& eps, doublereal& eps_dot,doublereal& delta,doublereal& SoD,doublereal& SoV,doublereal& my,doublereal& beta, doublereal r, doublereal alpha)const
{
        const Vec3& X1 = m_pShaft->GetXCurr();
        const Vec3& X2 = m_pBearing->GetXCurr();
        const Vec3& X1_dot = m_pShaft->GetVCurr();
        const Vec3& X2_dot = m_pBearing->GetVCurr();
        const Mat3x3& R1 = m_pShaft->GetRCurr();
        const Mat3x3& R2 = m_pBearing->GetRCurr();
        const Vec3& omega1 = m_pShaft->GetWCurr();
        const Vec3& omega2 = m_pBearing->GetWCurr();

		Vec3 o1_R1 = m_o1_R1;
    o1_R1(3) += r * m_bdat.b;
    Vec3 o2_R2 = m_o2_R2;
    o2_R2(3) += r * m_bdat.b;

    const Vec3 v_R2 = R2.MulTV( X1 - X2 + R1 * o1_R1 ) - o2_R2;
		const Vec3 d1_R2 = R2.MulTV(R1.GetCol(3));

		const doublereal lambda = m_lambda ? -v_R2(3) / d1_R2(3) : 0.;

		e_R2 = v_R2 + d1_R2 * lambda;
		const Vec3 v_dot_R2 = R2.MulTV( X1_dot - X2_dot + omega2.Cross(X2 - X1) + ( omega1 - omega2 ).Cross( R1 * o1_R1 ) );
		const Vec3 d1_dot_R2 = R2.MulTV( ( omega1 - omega2 ).Cross( R1.GetCol(3) ) );
		const doublereal lambda_dot = m_lambda ? -v_dot_R2(3) / d1_R2(3) + v_R2(3) / std::pow(d1_R2(3), 2) * d1_dot_R2(3) : 0.;

		// e_dot_R2 = R2^T * e_dot_I
    e_dot_R2 = R2.MulTV( X1_dot - X2_dot + omega1.Cross( R1 * o1_R1 ) - omega2.Cross( R2 * o2_R2 )
							+ R1.GetCol(3) * lambda_dot + omega1.Cross(R1.GetCol(3)) * lambda);
    const Vec3 l2_R2 = o2_R2 + e_R2;
		const Vec3 lambda_d1_R1 = Vec3(0.,0.,lambda);
		const Vec3 l1_I = R1 * ( o1_R1 + lambda_d1_R1 );

		omega_proj[0] = R2.GetCol(3).Dot(omega1);
		omega_proj[1] = R2.GetCol(3).Dot(omega2);

		const doublereal e_[2] = { e_R2(1), e_R2(2) }, e_dot_[2] = { e_dot_R2(1), e_dot_R2(2) };

		doublereal k[3];

    hydrodynamic_plain_bearing_force(m_bdat, omega_proj, e_, e_dot_, k, eps, eps_dot, delta, SoD, SoV, my, beta);

		F2_R2(1) = k[0];
		F2_R2(2) = k[1];
		F2_R2(3) = 0.;

		M2_R2(1) = 0.;
		M2_R2(2) = 0.;
		M2_R2(3) = k[2];

    alpha *= m_InitialAssemblyFactor.dGet();

		F2_R2 *= alpha;
		M2_R2 *= alpha;

		F2_I = R2 * F2_R2;
		M2_I = R2 * ( l2_R2.Cross( F2_R2 ) + M2_R2 );
		F1_I = -F2_I;
		M1_I = -l1_I.Cross( F2_I ) - R2 * M2_R2;
}


unsigned int
HydrodynamicPlainBearing::iGetNumPrivData(void) const
{
	return 0;
}

int
HydrodynamicPlainBearing::iGetNumConnectedNodes(void) const
{
	return 2; // 1x shaft + 1x bearing
}

void
HydrodynamicPlainBearing::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(iGetNumConnectedNodes());
	connectedNodes[0] = m_pShaft;
	connectedNodes[1] = m_pBearing;
}

void
HydrodynamicPlainBearing::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{

}

std::ostream&
HydrodynamicPlainBearing::Restart(std::ostream& out) const
{
	out << "hydrodynamic_plain_bearing_with_offset,\n"
			"		shaft," << m_pShaft->GetLabel() << ",\n"
			"		offset," << m_o1_R1 << ",\n"
			"		bearing," << m_pBearing->GetLabel() << ",\n"
			"		offset," << m_o2_R2 << "\n"
        "		bearing width," << m_bdat.b << ",\n"
        "		bearing diameter," << m_bdat.d << ",\n"
        "		relative clearance," << m_bdat.Psi << ",\n"
        "		oil viscosity," << m_bdat.eta << ",\n"
        "		initial assembly factor,";

    m_InitialAssemblyFactor.pGetDriveCaller()->Restart(out);

	out << ";\n";

	return out;
}

unsigned int
HydrodynamicPlainBearing::iGetInitialNumDof(void) const
{
	return 0;
}

void
HydrodynamicPlainBearing::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
HydrodynamicPlainBearing::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
HydrodynamicPlainBearing::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	WorkVec.ResizeReset(0);

	return WorkVec;
}

bool hydrodynamic_plain_bearing_set(void)
{
    UserDefinedElemRead *rf = new UDERead<HydrodynamicPlainBearing>;

	if (!SetUDE("hydrodynamic_plain_bearing_with_offset", rf))
	{
		delete rf;
		return false;
	}

	return true;
}

#ifndef STATIC_MODULES

extern "C"
{
int module_init(const char *module_name, void *pdm, void *php)
{
	if (!hydrodynamic_plain_bearing_set())
	{
		silent_cerr("hydrodynamic_plain_bearing: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

}

#endif // ! STATIC_MODULE

