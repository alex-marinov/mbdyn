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

// #define MATVEC_DEBUG 1
// #define GRADIENT_DEBUG 1

#include <dataman.h>
#include <userelem.h>
#include <matvecass.h>

#include "module-hydrodynamic_plain_bearing.h"

#if defined(USE_AUTODIFF)

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
    template <typename T>
    inline void
    AssRes(grad::GradientAssVec<T>& WorkVec,
           doublereal dCoef,
           const grad::GradientVectorHandler<T>& XCurr,
           const grad::GradientVectorHandler<T>& XPrimeCurr,
           enum grad::FunctionCall func);
    template <typename T>
    inline void
    InitialAssRes(grad::GradientAssVec<T>& WorkVec,
                  const grad::GradientVectorHandler<T>& XCurr,
                  enum grad::FunctionCall func);

	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
    virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
    virtual void Update(const VectorHandler& XCurr, 
                        const VectorHandler& XPrimeCurr);
	unsigned int iGetNumPrivData(void) const;
    virtual unsigned int iGetPrivDataIdx(const char *s) const;
    virtual doublereal dGetPrivData(unsigned int i) const;
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
    template <typename T>
    inline void
    UnivAssRes(grad::GradientAssVec<T>& WorkVec,
               doublereal dCoef,
               const grad::GradientVectorHandler<T>& XCurr,
               enum grad::FunctionCall func);
    
    template <typename T>
    struct OutputData {
        OutputData();
        
        doublereal r;
        doublereal alpha;
        grad::Vector<T, 3> e_R2;
        grad::Vector<T, 3> e_dot_R2;
        grad::Vector<T, 2> omega_proj;
        grad::Vector<T, 3> F2_R2;
        grad::Vector<T, 3> M2_R2;
        grad::Vector<T, 3> F2_I;
        grad::Vector<T, 3> M2_I;
        grad::Vector<T, 3> F1_I;
        grad::Vector<T, 3> M1_I;
        T eps;
        T eps_dot;
        T delta;
        T SoD;
        T SoV;
        T mu;
        T beta;
        bool bUpdated;
    };

    void UpdateOutput(OutputData<doublereal>& oOutput) const;
    
    template <typename T>
    inline
    void ComputeResidual(OutputData<T>& oOutput,
                         doublereal r,
                         doublereal alpha,
                         doublereal dCoef,
                         grad::FunctionCall eFunc,
                         grad::LocalDofMap* pDofMap) const;
    
    struct Bearing2D
    {
        doublereal b, D, Psi, eta, eps_max, s, a[9];

        Bearing2D();
        void Initialize();
        
        template <typename T>
        void SommerfeldNumbers(const T& eps,
                               const grad::Vector<T, 2>& omega,
                               const T& delta_dot, 
                               T& sod,
                               T& sov,
                               T& beta,
                               T& mu) const;

        template <typename T>
        void SommerfeldNumbersExt(const T& eps,
                                  const grad::Vector<T, 2>& omega,
                                  const T& delta_dot,
                                  T& sod,
                                  T& sov,
                                  T& beta, 
                                  T& mu) const;

        template <typename T>
        void UpdateBearingForce(OutputData<T>& oOutput) const;        
    };
    
    static const grad::index_type iNumADVars = 12;
	const StructNode* m_pShaft;
	const StructNode* m_pBearing;
	Vec3 m_o1_R1;
	Vec3 m_o2_R2;
    Bearing2D m_bdat;
	DriveOwner m_InitialAssemblyFactor;
    integer m_iNumGaussPoints;
	const doublereal* m_r;
	const doublereal* m_alpha;
    mutable OutputData<doublereal> m_output[6];
    
    integer m_iNumOutputPoints;
    bool m_lambda, m_bInitialAssembly;
    grad::LocalDofMap m_oDofMap;
	static const doublereal s_r1[1], s_alpha1[1];
	static const doublereal s_r2[2], s_alpha2[2];
	static const doublereal s_r3[3], s_alpha3[3];
	static const doublereal s_r6[6], s_alpha6[6];

    static const integer sm_iNumPrivData = 18;
    static const struct PrivData {
        char szName[13];
    } sm_rgPrivData[sm_iNumPrivData];
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

const HydrodynamicPlainBearing::PrivData HydrodynamicPlainBearing::sm_rgPrivData[sm_iNumPrivData] = {
    {"e1"},       // 0
    {"e2"},       // 1
    {"eP1"},      // 2
    {"eP2"},      // 3
    {"omega1"},   // 4
    {"omega2"},   // 5
    {"epsilon"},  // 6
    {"epsilonP"}, // 7
    {"delta"},    // 8
    {"F2x"},      // 9
    {"F2y"},      // 10
    {"M2z"},      // 11
    {"SoD"},      // 12
    {"SoV"},      // 13
    {"mu"},       // 14
    {"beta"},     // 15
    {"minh"},     // 16
    {"Pff"}       // 17
};

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
	m_lambda(true),
        m_bInitialAssembly(false)
{
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
            "		shaft, (label) <shaft_node>,\n"
			"		[offset, (Vec3) <o1>,]\n"
            "		bearing, (label) <bearing_node>,\n"
			"		[offset, (Vec3) <o2>,]\n"
            "		bearing width, (real) <b>,\n"
			"		{shaft diameter, (real) <d> | bearing_diameter, (real) <D>,}\n"
            "		relative clearance, (real) <Psi>,\n"
            "		oil viscosity, (real) <eta>,\n"
            "		initial assembly factor, (DriveCaller) <assembly_factor>,\n"
            "           [number of gauss points, (integer) <num_gauss_points>]\n"
            "           [output points, (integer) <num_output_points>\n"
            "               [, {gauss point, (integer) <index_1> | custom, r, (real) <r_1>, alpha, (real) <alpha_1>}]]\n"
            "       [extend shaft to bearing center, {yes | no | (bool) <extend>}]\n"
            "           [epsilon max, (real) <epsilon_max>]\n"
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

    doublereal d = -1., D = -1.;

	if (HP.IsKeyWord("shaft" "diameter"))
	{
        d = HP.GetReal();

        if (d <= 0)
        {
            silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): shaft diameter must be greater than zero at line " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
	}
	else
	{
		if ( !( HP.IsKeyWord("bearing_diameter") || HP.IsKeyWord("bearing" "diameter") ) )
		{
			silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"bearing diameter\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

        D = HP.GetReal();

        if (D <= 0)
        {
            silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): bearing diameter must be greater than zero at line " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
	}

	if ( !( HP.IsKeyWord("relative_clearance") || HP.IsKeyWord("relative" "clearance") ) )
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"relative clearance\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
        
    m_bdat.Psi = HP.GetReal();

    if (m_bdat.Psi <= 0 || m_bdat.Psi >= 1)
	{
        silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword relative clearance must be greater than zero and less than one at line " << HP.GetLineData() << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
    
    if (d > 0.)
    {
        m_bdat.D = d / (1 - m_bdat.Psi);
    }
    else
    {
        m_bdat.D = D;
	}

	if ( !( HP.IsKeyWord("oil_viscosity") || HP.IsKeyWord("oil" "viscosity") ) )
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"oil viscosity\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
        
    m_bdat.eta = HP.GetReal();

	m_InitialAssemblyFactor.Set(( HP.IsKeyWord("initial_assembly_factor") || HP.IsKeyWord("initial" "assembly" "factor") ) ? HP.GetDriveCaller() : new OneDriveCaller());

    m_iNumGaussPoints = HP.IsKeyWord("number" "of" "gauss" "points") ? HP.GetInt() : 1;

#define GAUSS_POINT_NUM_(num)                                      \
    num:                                                                \
        ASSERT(m_iNumGaussPoints == sizeof(s_r##num)/sizeof(s_r##num[0])); \
        ASSERT(m_iNumGaussPoints == sizeof(s_alpha##num)/sizeof(s_alpha##num[0])); \
        m_r = s_r##num;                                                 \
        m_alpha = s_alpha##num;

    switch (m_iNumGaussPoints)
		{
    case GAUSS_POINT_NUM_(1); break;
    case GAUSS_POINT_NUM_(2); break;
    case GAUSS_POINT_NUM_(3); break;
    case GAUSS_POINT_NUM_(6); break;
    default:
        silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel()
                    << "): integration rule with " << m_iNumGaussPoints
                    << " gauss points not supported at line "
                    << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

    m_bdat.Initialize();

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

    if (HP.IsKeyWord("initial" "assembly")) {
        m_bInitialAssembly = HP.GetYesNoOrBool();
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
        << m_bdat.D * (1. - m_bdat.Psi) << " "
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
            UpdateOutput(m_output[i]);

			// output the current state: the column layout is as follows
			// (all quantities are refered to the reference frame of the bearing node)

			// 1     2     3         4         5         6         7    8        9      10  11  12  13   14   15  16
            // e(1)  e(2)  e_dot(1)  e_dot(2)  omega(1)  omega(2)  eps  eps_dot  delta  Fx  Fy  Mz  SoD  SoV  mu  beta

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

            os 	<< m_output[i].e_R2(1) << ' ' 			// 1
                << m_output[i].e_R2(2) << ' ' 			// 2
                << m_output[i].e_dot_R2(1) << ' '		// 3
                << m_output[i].e_dot_R2(2) << ' '		// 4
                << m_output[i].omega_proj(1) << ' '		// 5
                << m_output[i].omega_proj(2) << ' '		// 6
                << m_output[i].eps << ' ' 			// 7
                << m_output[i].eps_dot << ' '			// 8
                << m_output[i].delta << ' ' 			// 9
                << m_output[i].F2_R2(1) << ' ' 			// 10
                << m_output[i].F2_R2(2) << ' ' 			// 11
                << m_output[i].M2_R2(3) << ' ' 			// 12
                << m_output[i].SoD << ' ' 			// 13
                << m_output[i].SoV << ' ' 			// 14
                << m_output[i].mu << ' ' 			// 15
                << m_output[i].beta << ' ';			// 16
		}

		os << std::endl;
	}
}

void
HydrodynamicPlainBearing::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
    *piNumRows = 12 * m_iNumGaussPoints;
	*piNumCols = 12;
}

template <typename T>
inline void
HydrodynamicPlainBearing::UnivAssRes(grad::GradientAssVec<T>& WorkVec,
	doublereal dCoef,
                                     const grad::GradientVectorHandler<T>& XCurr,
                                     enum grad::FunctionCall func)
{
    using namespace grad;

    typedef Matrix<T, 3, 3> Mat3x3;
    typedef Vector<T, 3> Vec3;

    const integer intShaftMomentumIndex = m_pShaft->iGetFirstMomentumIndex();
    const integer intBearingMomentumIndex = m_pBearing->iGetFirstMomentumIndex();

    OutputData<T> oOutput;

    for (integer i = 0; i < m_iNumGaussPoints; ++i)
	{
        ComputeResidual(oOutput,
                        m_r[i],
                        m_alpha[i],
                        dCoef,
                        func,
                        &m_oDofMap);

        // 1     2     3     4     5     6     7     8     9     10    11    12
        // F1(1) F1(2) F1(3) M1(1) M1(2) M1(3) F2(1) F2(2) F2(3) M2(1) M2(2) M2(3)
        WorkVec.AddItem(intShaftMomentumIndex + 1,  oOutput.F1_I);
        WorkVec.AddItem(intShaftMomentumIndex + 4,  oOutput.M1_I);
        WorkVec.AddItem(intBearingMomentumIndex + 1,  oOutput.F2_I);
        WorkVec.AddItem(intBearingMomentumIndex + 4, oOutput.M2_I);
    }
}


template <typename T>
inline void HydrodynamicPlainBearing::AssRes(grad::GradientAssVec<T>& WorkVec,
                                             doublereal dCoef,
                                             const grad::GradientVectorHandler<T>& XCurr,
                                             const grad::GradientVectorHandler<T>& XPrimeCurr,
                                             enum grad::FunctionCall func)
{
    UnivAssRes(WorkVec, dCoef, XCurr, func);
}

template <typename T>
inline void HydrodynamicPlainBearing::InitialAssRes(grad::GradientAssVec<T>& WorkVec,
                                                    const grad::GradientVectorHandler<T>& XCurr,
                                                    enum grad::FunctionCall func)
{
    UnivAssRes(WorkVec, 1., XCurr, func);
}


VariableSubMatrixHandler&
HydrodynamicPlainBearing::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
    using namespace grad;

    GradientAssVec<Gradient<iNumADVars> >::AssJac(this,
                                                  WorkMat.SetSparse(),
                                                  dCoef,
                                                  XCurr,
                                                  XPrimeCurr,
                                                  REGULAR_JAC,
                                                  &m_oDofMap);
    return WorkMat;
}


SubVectorHandler&
HydrodynamicPlainBearing::AssRes(SubVectorHandler& WorkVec,
                                 doublereal dCoef,
                                 const VectorHandler& XCurr,
                                 const VectorHandler& XPrimeCurr)
{
    using namespace grad;

    GradientAssVec<doublereal>::AssRes(this,
                                       WorkVec,
                                       dCoef,
                                       XCurr,
                                       XPrimeCurr,
                                       grad::REGULAR_RES);

    return WorkVec;
}

void HydrodynamicPlainBearing::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
    Update(X, XP);
}

void HydrodynamicPlainBearing::Update(const VectorHandler& XCurr, 
                                      const VectorHandler& XPrimeCurr)
{
    for (integer i = 0; i < m_iNumOutputPoints; ++i) {
        m_output[i].bUpdated = false;
    }
}

void HydrodynamicPlainBearing::UpdateOutput(OutputData<doublereal>& oOutput) const
{
    if (oOutput.bUpdated) {
        return;
    }

    ComputeResidual(oOutput, oOutput.r, oOutput.alpha, 1., grad::REGULAR_RES, nullptr);

    oOutput.bUpdated = true;
}

template <typename T>
void HydrodynamicPlainBearing::ComputeResidual(OutputData<T>& oOutput,
                                               doublereal r,
                                               doublereal alpha,
                                               doublereal dCoef,
                                               grad::FunctionCall eFunc,
                                               grad::LocalDofMap* pDofMap) const
{
    using namespace grad;
    typedef Vector<T, 3> Vec3;
    typedef Matrix<T, 3, 3> Mat3x3;

    Vec3 X1, X2, X1_dot, X2_dot, omega1, omega2;
    Mat3x3 R1, R2;

    m_pShaft->GetXCurr(X1, dCoef, eFunc, pDofMap);
    m_pBearing->GetXCurr(X2, dCoef, eFunc, pDofMap);
    m_pShaft->GetVCurr(X1_dot, dCoef, eFunc, pDofMap);
    m_pBearing->GetVCurr(X2_dot, dCoef, eFunc, pDofMap);
    m_pShaft->GetRCurr(R1, dCoef, eFunc, pDofMap);
    m_pBearing->GetRCurr(R2, dCoef, eFunc, pDofMap);

    m_pShaft->GetWCurr(omega1, dCoef, eFunc, pDofMap);
    m_pBearing->GetWCurr(omega2, dCoef, eFunc, pDofMap);

    Vec3 o1_R1(m_o1_R1);

    o1_R1(3) += r * m_bdat.b;

    Vec3 o2_R2(m_o2_R2);

    o2_R2(3) += r * m_bdat.b;

    const Vec3 v_R2 = Transpose(R2) * Vec3(X1 - X2 + R1 * o1_R1) - o2_R2;
    const Vec3 d1_R2 = Transpose(R2) * R1.GetCol(3);

    T lambda(0.);

    if (m_lambda) {
        lambda = -v_R2(3) / d1_R2(3);
    }

    oOutput.e_R2 = v_R2 + d1_R2 * lambda;

    const Vec3 v_dot_R2 = Transpose(R2) * Vec3(X1_dot - X2_dot + Cross(omega2, X2 - X1) + Cross(omega1 - omega2, R1 * o1_R1));
    const Vec3 d1_dot_R2 = Transpose(R2) * Vec3(Cross(omega1 - omega2, R1.GetCol(3)));
    T lambda_dot(0.);

    if (m_lambda) {
        lambda_dot = -v_dot_R2(3) / d1_R2(3) + v_R2(3) / pow(d1_R2(3), 2) * d1_dot_R2(3);
    }

    // e_dot_R2 = R2^T * e_dot_I
    oOutput.e_dot_R2 = Transpose(R2) * Vec3( X1_dot - X2_dot + Cross(omega1,  R1 * o1_R1 ) - Cross(omega2, R2 * o2_R2)
                                     + R1.GetCol(3) * lambda_dot + Cross(omega1, R1.GetCol(3)) * lambda);
    const Vec3 l2_R2 = o2_R2 + oOutput.e_R2;
    const Vec3 lambda_d1_R1(T(0.), T(0.), lambda);
    const Vec3 l1_I = R1 * Vec3( o1_R1 + lambda_d1_R1 );

    oOutput.omega_proj(1) = Dot(R2.GetCol(3), omega1);
    oOutput.omega_proj(2) = Dot(R2.GetCol(3), omega2);

    if (typeid(T) != typeid(doublereal)) {
        if (oOutput.e_R2(1) == 0. && oOutput.e_R2(2) == 0.) {
            oOutput.e_R2(1) += m_bdat.s * std::numeric_limits<doublereal>::epsilon();
        }

        if (oOutput.e_dot_R2(1) == 0. && oOutput.e_dot_R2(2) == 0.) {
            oOutput.e_dot_R2(1) += m_bdat.s * std::numeric_limits<doublereal>::epsilon();
        }
    }

    m_bdat.UpdateBearingForce(oOutput);

    alpha *= m_InitialAssemblyFactor.dGet();

    oOutput.F2_R2 *= alpha;
    oOutput.M2_R2 *= alpha;

    oOutput.F2_I = R2 * oOutput.F2_R2;
    oOutput.M2_I = R2 * Vec3(Cross(l2_R2, oOutput.F2_R2) + oOutput.M2_R2);
    oOutput.F1_I = -oOutput.F2_I;
    oOutput.M1_I = -Cross(l1_I, oOutput.F2_I) - R2 * oOutput.M2_R2;
}


unsigned int
HydrodynamicPlainBearing::iGetNumPrivData(void) const
{
    return sm_iNumPrivData * m_iNumOutputPoints;
}

unsigned int HydrodynamicPlainBearing::iGetPrivDataIdx(const char *s) const
{
    const char* pszIndex = strchr(s, '[');

    if (!pszIndex) {
        return 0;
    }

    int iPrivData;

    for (iPrivData = 0; iPrivData < sm_iNumPrivData; ++iPrivData) {
        if (0 == strncmp(s, sm_rgPrivData[iPrivData].szName, pszIndex - s)) {
            break;
        }
    }

    if (iPrivData >= sm_iNumPrivData) {
        return 0;
	}

    int iOutputLoc;

    if (1 != sscanf(++pszIndex, "%d", &iOutputLoc)) {
        return 0;
    }

    if (iOutputLoc < 1 || iOutputLoc > m_iNumOutputPoints) {
        return 0;
	}

    return sm_iNumPrivData * (iOutputLoc - 1) + iPrivData + 1;
}

doublereal HydrodynamicPlainBearing::dGetPrivData(unsigned int i) const
{
    div_t oIndex = div(i - 1, sm_iNumPrivData);

    if (oIndex.quot < 0 || oIndex.quot >= m_iNumOutputPoints) {
        silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): invalid private data index " << i << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

    UpdateOutput(m_output[oIndex.quot]);

    switch (oIndex.rem) {
    case 0:
    case 1:
        return m_output[oIndex.quot].e_R2(oIndex.rem + 1);
    case 2:
    case 3:
        return m_output[oIndex.quot].e_dot_R2(oIndex.rem - 1);
    case 4:
    case 5:
        return m_output[oIndex.quot].omega_proj(oIndex.rem - 3);
    case 6:
        return m_output[oIndex.quot].eps;
    case 7:
        return m_output[oIndex.quot].eps_dot;
    case 8:
        return m_output[oIndex.quot].delta;
    case 9:
    case 10:
        return m_output[oIndex.quot].F2_R2(oIndex.rem - 8);
    case 11:
        return m_output[oIndex.quot].M2_R2(3);
    case 12:
        return m_output[oIndex.quot].SoD;
    case 13:
        return m_output[oIndex.quot].SoV;
    case 14:
        return m_output[oIndex.quot].mu;
    case 15:
        return m_output[oIndex.quot].beta;
    case 16:
        return 0.5 * (1. - fabs(m_output[oIndex.quot].eps)) * m_bdat.s;
    case 17:
        return (m_output[oIndex.quot].omega_proj(2) - m_output[oIndex.quot].omega_proj(1)) * m_output[oIndex.quot].M2_R2(3);
    default:
        silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): invalid private data index " << i << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
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
        "		bearing diameter," << m_bdat.D << ",\n"
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
    return 0u;
}

void
HydrodynamicPlainBearing::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
    if (m_bInitialAssembly) {
        *piNumRows = 12 * m_iNumGaussPoints;
        *piNumCols = 24;
    } else {
        *piNumRows = *piNumCols = 0;
    }
}

VariableSubMatrixHandler&
HydrodynamicPlainBearing::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
    using namespace grad;

    if (m_bInitialAssembly) {
        GradientAssVec<Gradient<2 * iNumADVars> >::InitialAssJac(this,
                                                                 WorkMat.SetSparse(),
                                                                 XCurr,
                                                                 grad::INITIAL_ASS_JAC,
                                                                 &m_oDofMap);
    } else {
	WorkMat.SetNullMatrix();
    }

	return WorkMat;
}

SubVectorHandler&
HydrodynamicPlainBearing::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
    using namespace grad;

    if (m_bInitialAssembly) {
        GradientAssVec<doublereal>::InitialAssRes(this,
                                                  WorkVec,
                                                  XCurr,
                                                  grad::INITIAL_ASS_RES);
    } else {
	WorkVec.ResizeReset(0);
    }

	return WorkVec;
}

template <typename T>
HydrodynamicPlainBearing::OutputData<T>::OutputData()
    :r(0.),
     alpha(0.),
     eps(0.),
     eps_dot(0.),
     delta(0.),
     SoD(0.),
     SoV(0.),
     mu(0.),
     beta(0.),
     bUpdated(false)
{
    
}

/* hydrodynamic_plain_bearing_force.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/* Subroutine */
HydrodynamicPlainBearing::Bearing2D::Bearing2D()
{
    std::memset(this, 0, sizeof(*this));
}

void HydrodynamicPlainBearing::Bearing2D::Initialize()
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Function Body */
    s = D * Psi;
/*     coefficients for SoD for rotation according to Butenschoen */
/* Computing 2nd power */
    d__1 = b / D;
/* Computing 3rd power */
    d__2 = b / D;
/* Computing 4th power */
    d__3 = d__2 * d__2;

    a[1 - 1] = 1.1642 - b / D * 1.9456 + d__1 * d__1 * 7.1161 - d__2 * (d__2 *
                                                                        d__2) * 10.1073 + d__3 * d__3 * 5.0141;    
    a[2 - 1] = -1.000026 - b / D * .023634 - d__1 * d__1 * .4215 - d__2 * (
        d__2 * d__2) * .038817 - d__3 * d__3 * .090551;
/*     coefficients for beta for rotation */
    a[3 - 1] = 1.152624 - b / D * .104565;
    a[4 - 1] = b / D * .798745 - 2.5905;
    a[5 - 1] = 8.73393 - b / D * 2.3291;
    a[6 - 1] = b / D * 3.424337 - 13.3414;
    a[7 - 1] = 6.6294 - b / D * 1.591732;
/*     coefficients for displacement */
    a[8 - 1] = b / D * 3.2415 + .70038 - d__1 * d__1 * 12.2486 + d__2 * (d__2 
                                                                         * d__2) * 18.895 - d__3 * d__3 * 9.3561;
    a[9 - 1] = b / D * .0157434 - .999935 - d__1 * d__1 * .74224 + d__2 * (
        d__2 * d__2) * .42278 - d__3 * d__3 * .368928;
} /* hydrodynaic_plain_bearing_init__ */

/* Subroutine */
template <typename T>
void HydrodynamicPlainBearing::Bearing2D::SommerfeldNumbers(const T& eps,
                                                            const grad::Vector<T, 2>& omega,
                                                            const T& delta_dot, 
                                                            T& sod,
                                                            T& sov,
                                                            T& beta,
                                                            T& mu) const
{
    static const doublereal c_b3 = -2.5;
    /* System generated locals */
    T d__1, d__2, d__3, d__4, d__5;
/*     Sommerfeld number for rotation according to Butenschoen 1976 */
    /* Parameter adjustments */
    /* Function Body */
/* Computing 2nd power */
    d__1 = b / D;
/* Computing 2nd power */
    d__3 = eps;
/* Computing 2nd power */
    d__2 = 1. - d__3 * d__3;
/* Computing 2nd power */
    d__4 = eps;
/* Computing 2nd power */
    d__5 = eps;
    sod = d__1 * d__1 * fabs(eps) / (d__2 * d__2 * 2.) * sqrt((1. - d__4 * 
	    d__4) * 9.8696044010893385 + d__5 * d__5 * 16.) * a[1 - 1] * (fabs(
	    eps) - 1.) / (a[2 - 1] + fabs(eps));
/*     Sommerfeld number for displacement according to Butenschoen 1976 */
/* Computing 2nd power */
    d__1 = b / D;
/* Computing 2nd power */
    d__3 = eps;
    d__2 = 1. - d__3 * d__3;
/* Computing 2nd power */
    d__4 = eps;
/* Computing 2nd power */
    d__5 = eps;
    sov = d__1 * d__1 * 4. * pow(d__2, c_b3) * ((M_PI/2. - 
	    acos(eps) * .5) * (d__4 * d__4 * 2. + 1.) + eps * 1.5 * sqrt(1. 
	    - d__5 * d__5)) * a[8 - 1] * (1. - eps) / (-a[9 - 1] - eps);
/*     angle between force for rotation and minimum clearance according to Butenschoen 1976 */
/* Computing 2nd power */
    d__1 = eps;
/* Computing 2nd power */
    d__2 = eps;
/* Computing 3rd power */
    d__3 = fabs(eps);
/* Computing 4th power */
    d__4 = eps * eps;
    
    beta = atan2(sqrt(1. - d__1 * d__1) * M_PI, fabs(eps) * 2.) *
	     (a[3 - 1] + a[4 - 1] * fabs(eps) + a[5 - 1] * (d__2 * d__2) + a[6 - 1] * (d__3 * 
	    (d__3 * d__3)) + a[7 - 1] * (d__4 * d__4));
    
    if (fabs(eps) < 1e-6) {
/*     avoid division infinite by infinite in case of zero relative eccentricity */
/*     use analytical limit of abs_MR for eps going to zero */
	mu = std::numeric_limits<doublereal>::max();
    } else {
/*     friction coefficient according to Butenschoen */
/* Computing 2nd power */
	mu = Psi * (fabs((omega(1) - omega(2)) / (omega(2) + omega(1) - 
		delta_dot * 2.)) * M_PI / (sqrt(1. - 
		eps * eps) * sod) + sin(beta) * fabs(eps) / 2.);
    }
} /* SommerfeldNumbers */

/* Subroutine */ 
template <typename T>
void HydrodynamicPlainBearing::Bearing2D::SommerfeldNumbersExt(const T& eps,
                                                               const grad::Vector<T, 2>& omega,
                                                               const T& delta_dot,
                                                               T& sod,
                                                               T& sov,
                                                               T& beta, 
                                                               T& mu) const
{
    /* Function Body */
    if (fabs(eps) < eps_max) {
/*     According to the thesis of Butenschoen those approximations are valid until eps = 0.999 */
	SommerfeldNumbers(eps,
                          omega, 
                          delta_dot,
                          sod,
                          sov,
                          beta,
                          mu);
    } else {
        /* Local variables */        
        T sod1, eps0, eps1, sov1, beta1, mu1;
/*     Do a linear extrapolation above eps_max */
	eps0 = copysign(eps_max, eps);
        
	SommerfeldNumbers(eps0,
                          omega, 
                          delta_dot,
                          sod,
                          sov,
                          beta,
                          mu);
        
	eps1 = eps0 * (1. - sqrt(std::numeric_limits<doublereal>::epsilon()));
        
	SommerfeldNumbers(eps1,
                          omega, 
                          delta_dot,
                          sod1,
                          sov1,
                          beta1,
                          mu1);
        
	sod += T((sod1 - sod) / (eps1 - eps0) * (eps - eps0));
	sov += T((sov1 - sov) / (eps1 - eps0) * (eps - eps0));
    }
} /* SommerfeldNumbersExt */

/* Subroutine */ 
template <typename T>
void HydrodynamicPlainBearing::Bearing2D::UpdateBearingForce(OutputData<T>& oOutput) const
{
    static const doublereal c_b6 = 1.;
    /* System generated locals */
    T d__1, d__2;

    /* Local variables */
    T abs_e_dot__, delta_dot__, omega_res__, phi, abs_e__, 
        alpha, kappa, abs_fd__, abs_fv__, abs_mr__;
/* ----------------------------------------------------------------------------------------------------------
--------------- */
/*     hydrodynamic plain bearing calculation according to Butenschoen's theory */
/* ----------------------------------------------------------------------------------------------------------
--------------- */
/*     COORDINATE SYSTEM: */
/*     x ... axial direction */
/*     y, z ... radial direction */
/* ----------------------------------------------------------------------------------------------------------
--------------- */
/*     INPUT PARAMETERS */
/* ----------------------------------------------------------------------------------------------------------
--------------- */
/*     b   ... bearing width [m] */
/*     d   ... shaft diameter [m] */
/*     D   ... bearing diameter [m] */
/*     Psi ... relative radial clearance Psi = ( D - d ) / D [1] */
/*     eta ... dynamic oil viscosity [Pa*s] */
/*     omega(1) ... angular velocity of the shaft [rad/s] */
/*     omega(2) ... angular velocity of the bearing [rad/s] */
/*     e ...  radial eccentricity of the shaft */
/*     e(1) = ey  [m] */
/*     e(2) = ez  [m] */
/*     e_dot(1) ... velocity of the shaft relative to the bearing */
/*     e_dot(1) = ey_dot [m/s] */
/*     e_dot(2) = ez_dot [m/s] */

/* ----------------------------------------------------------------------------------------------------------
---------------- */
/*     OUTPUT PARAMETERS */
/* ----------------------------------------------------------------------------------------------------------
---------------- */
/*     F2_R2 ... force on the bearing */
/*     M2_R2 ... torque at the bearing */
/*     F2_R2(1) = Fy [N] */
/*     F2_R2(2) = Fz [N] */
/*     M2_R2(3) = Mx [Nm] */
/*     eps ...             relative eccentricity of the shaft [1] */
/*     eps_dot ...      time derivative of the relative eccentricity [1/s] */
/*     delta ...          angle of minimum clearance between shaft and bearing [rad] */
/*     SoD ....           Sommerfeld number for rotation [1] */
/*     SoV ...            Sommerfeld number for displacement [1] */
/*     mu ...             friction coefficient [N/N] */
/*     beta  ...          angle between force for rotation and minimum clearance [rad] */
/*     angle of the position with minimum clearance between shaft and bearing */
    
    /* Function Body */
    oOutput.delta = atan2(oOutput.e_R2(2), oOutput.e_R2(1));
/*     angle of the velocity vector of the shaft relative to the bearing */
    phi = atan2(oOutput.e_dot_R2(2), oOutput.e_dot_R2(1));
/*     angle between velocity vector and minimum clearance */
    kappa = phi - oOutput.delta;
/*     absolute value of the eccentricity of the shaft inside the bearing */
/* Computing 2nd power */
    d__1 = oOutput.e_R2(1);
/* Computing 2nd power */
    d__2 = oOutput.e_R2(2);
    abs_e__ = sqrt(d__1 * d__1 + d__2 * d__2);
/*     absolute value of the velocity of the shaft relative to the bearing */
/* Computing 2nd power */
    d__1 = oOutput.e_dot_R2(1);
/* Computing 2nd power */
    d__2 = oOutput.e_dot_R2(2);
    abs_e_dot__ = sqrt(d__1 * d__1 + d__2 * d__2);
/*     time derivative of the relative eccentricity of the shaft */
    oOutput.eps_dot = cos(kappa) * 2. * abs_e_dot__ / s;
/*     relative eccentricity of the shaft */
    oOutput.eps = abs_e__ * 2. / s;
    
    if (oOutput.eps_dot != 0.) {
/*     eps is positive if it's time derivative is positive too */
/*     attention the signum function is zero if eps_dot is zero */
/*     but eps must not be zero in this case */
	oOutput.eps *= copysign(c_b6, oOutput.eps_dot);
    }
/*     time derivative of angle of minimum clearance */
    if (abs_e__ == 0.) {
/*     avoid division by zero */
	delta_dot__ = 0.;
    } else {
/* Computing 2nd power */
	d__1 = oOutput.e_R2(2);
/* Computing 2nd power */
	d__2 = oOutput.e_R2(1);
	delta_dot__ = (oOutput.e_R2(1) * oOutput.e_dot_R2(2) - oOutput.e_R2(2) * oOutput.e_dot_R2(1)) / (d__1 * d__1 
		+ d__2 * d__2);
    }
    
    SommerfeldNumbersExt(oOutput.eps,
                         oOutput.omega_proj,
                         delta_dot__,
                         oOutput.SoD,
                         oOutput.SoV,
                         oOutput.beta,
                         oOutput.mu);
    
/*     effective hydrodynamic angular velocity according to Butenschoen 1976 */
    omega_res__ = oOutput.omega_proj(1) + oOutput.omega_proj(2) - delta_dot__ * 2.;
/*     angle of the force for rotation */
    alpha = oOutput.delta - oOutput.beta * copysign(c_b6, omega_res__);
/*     absolute value of the force for rotation */
/* Computing 2nd power */
    d__1 = Psi;
    abs_fd__ = oOutput.SoD * (b * D * eta * fabs(omega_res__)) / (d__1 * d__1);
/*     absolute value of the force for displacement */
/* Computing 2nd power */
    d__1 = Psi;
    abs_fv__ = oOutput.SoV * (b * D * eta * oOutput.eps_dot) / (d__1 * d__1);
    
    if (oOutput.mu >= std::numeric_limits<doublereal>::max()) {
/* Computing 2nd power */
	abs_mr__ = b * M_PI * (D * D) * eta * (oOutput.omega_proj(1) - oOutput.omega_proj(2)) / Psi / 2.;
    } else {
	abs_mr__ = oOutput.mu * abs_fd__ * D / 2. * copysign(c_b6, oOutput.omega_proj(1) - oOutput.omega_proj(2));
    }

/*     sum of force for rotation and force for displacement */
    oOutput.F2_R2(1) = abs_fd__ * cos(alpha) + abs_fv__ * cos(oOutput.delta);
    oOutput.F2_R2(2) = abs_fd__ * sin(alpha) + abs_fv__ * sin(oOutput.delta);
    oOutput.F2_R2(3) = 0.;
    
    oOutput.M2_R2(1) = 0.;
    oOutput.M2_R2(2) = 0.;
    oOutput.M2_R2(3) = abs_mr__;

} /* BearingForce */

#endif

bool hydrodynamic_plain_bearing_set(void)
{
#if defined(USE_AUTODIFF)
    UserDefinedElemRead *rf = new UDERead<HydrodynamicPlainBearing>;

	if (!SetUDE("hydrodynamic_plain_bearing_with_offset", rf))
	{
		delete rf;
		return false;
	}

	return true;
#else
    return false;
#endif
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

