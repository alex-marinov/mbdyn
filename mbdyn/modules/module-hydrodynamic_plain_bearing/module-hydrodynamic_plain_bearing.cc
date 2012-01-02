/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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
 AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
        Copyright (C) 2011(-2012) all rights reserved.

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

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>
#include <userelem.h>

#include "module-hydrodynamic_plain_bearing.h"
#include "hydrodynamic_plain_bearing_force.h"

class hydrodynamic_plain_bearing_with_offset: virtual public Elem, public UserDefinedElem
{
public:
	hydrodynamic_plain_bearing_with_offset(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~hydrodynamic_plain_bearing_with_offset(void);
	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		   doublereal dCoef,
		   const VectorHandler& XCurr,
		   const VectorHandler& XPrimeCurr);
	void ComputeResidual(Vec3& e_R2,
						 Vec3& e_dot_R2,
						 double omega_proj[2],
						 Vec3& F2_R2,
						 Vec3& M2_R2,
						 Vec3& F2_I,
						 Vec3& M2_I,
						 Vec3& F1_I,
						 Vec3& M1_I,
						 double& eps,
						 double& eps_dot,
						 double& delta,
						 double& SoD,
						 double& SoV,
						 double& my,
						 double& beta)const;
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
	StructNode* m_pShaft;
	StructNode* m_pBearing;
	Vec3 m_o1_R1;
	Vec3 m_o2_R2;

	doublereal m_b;		// bearing width
	doublereal m_d;		// shaft diameter
	doublereal m_Psi;	// relative clearance Psi = s / d; s = D - d
	doublereal m_eta;	// dynamic oil viscosity
	DriveOwner m_InitialAssemblyFactor;
	bool m_fcontact;	// use the contact model
	doublereal m_sP;	// contact stiffness when epsilon == 1
	doublereal m_DL;	// Lehr damping coefficient when epsilon == 1
	doublereal m_m;		// equivalent dynamic mass for which m_DL is valid for
	doublereal m_abs_FPrs1;	// contact force for epsilon == 1
	doublereal m_myP;	// coulomb friction coefficient
	doublereal m_signum_delta_omega; // maximum angular velocity difference where the friction torque is zero
};

hydrodynamic_plain_bearing_with_offset::hydrodynamic_plain_bearing_with_offset(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: 	Elem(uLabel, flag(0)),
	UserDefinedElem(uLabel, pDO),
	m_pShaft(0),
	m_pBearing(0),
	m_o1_R1(0.,0.,0.),
	m_o2_R2(0.,0.,0.),
	m_b(0.),
	m_d(0.),
	m_Psi(0.),
	m_eta(0.),
	m_InitialAssemblyFactor(0),
	m_fcontact(false),
	m_sP(0.),
	m_DL(0.),
	m_m(1.),
	m_abs_FPrs1(0.),
	m_myP(0.),
	m_signum_delta_omega(0.)
{
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
			"		bearing_diameter, (real) <d>,\n"
			"		relative_clearance, (real) <Psi>,\n"
			"		oil_viscosity, (real) <eta>,\n"
			"		initial_assembly_factor, (DriveCaller),\n"
			"           [contact,\n"
			"			stiffness, (real) <sP>,\n"
			"			damping_ratio, (real) <DL>,\n"
			"			rated_mass, (real) <m>,\n"
			"			force, (real) <abs_FPrs1>,\n"
			"			coulomb_friction_coefficient, (real) <myP>,\n"
			"			angular_velocity_threshold, (real) <signum_delta_omega>]\n"
			"\n"
			"   b ... bearing width [m]\n"
			"	d ... bearing diameter [m]\n"
			"	Psi ... relative clearance Psi = ( D - d ) / d [m/m]\n"
			"	eta ... dynamic oil viscosity [Pa*s]\n"
			"\n"
			"	sP ... contact stiffness for wall contact at epsilon == 1 [N/m]\n"
			"	DL ... Lehr damping ratio for wall contact at epsilon == 1 [1]\n"
			"	m ... equivalent dynamic mass used for damping force calculation  [kg]\n"
			"	abs_FPrs1 ... contact force at epsilon == 1	[N]\n"
			"	myP ... coulomb friction coefficient [N/N]\n"
			"	signum_delta_omega ... smallest angular velocity where coulomb friction is considered [rad/s]\n"
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
	m_pShaft = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	ASSERT(m_pShaft != 0);

	if ( HP.IsKeyWord("offset") )
		m_o1_R1 = HP.GetPosRel(ReferenceFrame(m_pShaft));

	if ( !HP.IsKeyWord("bearing"))
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"bearing\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_pBearing = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	ASSERT( m_pBearing != 0 );

	if ( HP.IsKeyWord("offset") )
		m_o2_R2 = HP.GetPosRel(ReferenceFrame(m_pBearing));

	if ( !( HP.IsKeyWord("bearing_width") || HP.IsKeyWord("bearing" "width") ) )
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"bearing width\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}		
	m_b = HP.GetReal();

	if ( !( HP.IsKeyWord("bearing_diameter") || HP.IsKeyWord("bearing" "diameter") ) )
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"bearing diameter\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_d = HP.GetReal();

	if ( !( HP.IsKeyWord("relative_clearance") || HP.IsKeyWord("relative" "clearance") ) )
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"relative clearance\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_Psi = HP.GetReal();

	if ( !( HP.IsKeyWord("oil_viscosity") || HP.IsKeyWord("oil" "viscosity") ) )
	{
		silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"oil viscosity\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_eta = HP.GetReal();

	m_InitialAssemblyFactor.Set(( HP.IsKeyWord("initial_assembly_factor") || HP.IsKeyWord("initial" "assembly" "factor") ) ? HP.GetDriveCaller() : new OneDriveCaller());

	if ( HP.IsKeyWord("contact") )
	{
		m_fcontact = true;

		if ( !HP.IsKeyWord("stiffness") )
		{
			silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"stiffness\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_sP = HP.GetReal();

		if ( !( HP.IsKeyWord("damping_ratio") || HP.IsKeyWord("damping" "ratio") ) )
		{
			silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"damping ratio\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_DL = HP.GetReal();

		if ( !( HP.IsKeyWord("rated_mass") || HP.IsKeyWord("rated" "mass") ) )
		{
			silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"rated mass\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_m = HP.GetReal();

		if ( !HP.IsKeyWord("force") )
		{
			silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"force\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_abs_FPrs1 = HP.GetReal();

		if ( !( HP.IsKeyWord("coulomb_friction_coefficient") || HP.IsKeyWord("coulomb" "friction" "coefficient") ) )
		{
			silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"coulomb friction coefficient\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_myP = HP.GetReal();

		if ( !( HP.IsKeyWord("angular velocity threshold") || HP.IsKeyWord("angular" "velocity" "threshold") ) )
		{
			silent_cerr("hydrodynamic_plain_bearing_with_offset(" << GetLabel() << "): keyword \"angular velocity threshold\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_signum_delta_omega = HP.GetReal();
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	std::ostream& out = pDM->GetLogFile();

	out << "hydrodynamic_plain_bearing_with_offset: "
		<< uLabel << " "
		<< m_pShaft->GetLabel() << " "
		<< m_o1_R1 << " "
		<< m_pBearing->GetLabel() << " "
		<< m_o2_R2 << " "
		<< m_b << " "
		<< m_d << " "
		<< m_Psi << " "
		<< m_eta << " ";

	if ( m_fcontact )
	{
		out << m_sP << " "
			<< m_DL << " "
			<< m_m << " "
			<< m_abs_FPrs1 << " "
			<< m_myP << " "
			<< m_signum_delta_omega << " ";
	}

	out << std::endl;
}

hydrodynamic_plain_bearing_with_offset::~hydrodynamic_plain_bearing_with_offset(void)
{
	// destroy private data
	#ifdef DEBUG
		std::cerr << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << std::endl;
	#endif
}

void
hydrodynamic_plain_bearing_with_offset::Output(OutputHandler& OH) const
{
	if ( fToBeOutput() )
	{
		double eps, eps_dot, delta, SoD, SoV, my, beta;

		Vec3 e_R2, e_dot_R2;
		double omega_proj[2];
		Vec3 F2_R2, M2_R2;
		Vec3 F2_I, M2_I, F1_I, M1_I;

		ComputeResidual(e_R2,e_dot_R2,omega_proj,F2_R2,M2_R2,F2_I,M2_I,F1_I,M1_I,eps,eps_dot,delta,SoD,SoV,my,beta);

		if ( OH.UseText(OutputHandler::LOADABLE) )
		{
			std::ostream& os = OH.Loadable();

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

			os 	<< std::setw(8)
				<< GetLabel() << ' ' 		// 0
				<< e_R2(1) << ' ' 			// 1
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
				<< beta << ' '				// 16
				<< std::endl;
		}
	}
}

void
hydrodynamic_plain_bearing_with_offset::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 12;
}

VariableSubMatrixHandler&
hydrodynamic_plain_bearing_with_offset::AssJac(VariableSubMatrixHandler& WorkMatV,
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

	const Vec3 v_R2 = R2.MulTV( X1 - X2 + R1 * m_o1_R1 ) - m_o2_R2;
	const Vec3 d1_R2 = R2.MulTV(R1.GetCol(3));
	const doublereal lambda = -v_R2(3) / d1_R2(3);
	Vec3 e_R2 = v_R2 + d1_R2 * lambda;

	// FIXME: nan values if e_R2(1) == 0 && e_R2(2) == 0
	if ( e_R2(1) == 0.0 && e_R2(2) == 0.0 )
			e_R2(1) = std::numeric_limits<doublereal>::epsilon() * m_d * m_Psi / 2.;

	const Vec3 v_dot_R2 = R2.MulTV( X1_dot - X2_dot + omega2.Cross(X2 - X1) + ( omega1 - omega2 ).Cross( R1 * m_o1_R1 ) );
	const Vec3 d1_dot_R2 = R2.MulTV( ( omega1 - omega2 ).Cross( R1.GetCol(3) ) );
	const doublereal lambda_dot = -v_dot_R2(3) / d1_R2(3) + v_R2(3) / pow(d1_R2(3),2) * d1_dot_R2(3);

	Vec3 e_dot_R2 = R2.MulTV( X1_dot - X2_dot + omega1.Cross( R1 * m_o1_R1 ) - omega2.Cross( R2 * m_o2_R2 )
							+ R1.GetCol(3) * lambda_dot + omega1.Cross(R1.GetCol(3)) * lambda);

	// FIXME: nan values if e_dot_R2(1) == 0 && e_dot_R2(2) == 0
	if ( e_dot_R2(1) == 0.0 && e_dot_R2(2) == 0.0 )
			e_dot_R2(1) = std::numeric_limits<doublereal>::epsilon() * m_d * m_Psi / 2.;

	const Vec3 lambda_d1_R1 = Vec3(0.,0.,lambda);
	const Vec3 l1_R1 = m_o1_R1 + lambda_d1_R1;
	const Vec3 l1_I = R1 * l1_R1;
	const Vec3 l2_R2 = m_o2_R2 + e_R2;
	const double omega_proj[2] = { R2.GetCol(3).Dot(omega1),
				       R2.GetCol(3).Dot(omega2) };

	const double e_[2] = { e_R2(1), e_R2(2) }, e_dot_[2] = { e_dot_R2(1), e_dot_R2(2) };

	//												  |	0   1   2   3   4   5  |
	static const double ed[2][NBDIRSMAX] 	 	  = { { 1., 0., 0., 0., 0., 0. },	// 0 | inner derivative of the eccentricity of the shaft in direction 1 of the reference frame of the bearing (R2)
						     	 	 	 	 	 	  { 0., 1., 0., 0., 0., 0. } };	// 1 | inner derivative of the eccentricity of the shaft in direction 2 of the reference frame of the bearing (R2)
	static const double e_dotd[2][NBDIRSMAX] 	  = { { 0., 0., 1., 0., 0., 0. },	// 0 | inner derivative of relative velocity of the shaft in direction 1 of the reference frame of the bearing (R2)
						     	 	 	 	 	 	  { 0., 0., 0., 1., 0., 0. } }; // 1 | inner derivative of relative velocity of the shaft in direction 2 of the reference frame of the bearing (R2)
	static const double omega_projd[2][NBDIRSMAX] = { { 0., 0., 0., 0., 1., 0. },	// 0 | inner derivative of the angular velocity of the shaft in direction 3 of the reference frame of the bearing (R2)
												 	  { 0., 0., 0., 0., 0., 1. } }; // 1 | inner derivative of the angular velocity of the bearing in direction 3 of the reference frame of the bearing (R2)

	double k[3];			// force vector at the bearing (not used for the evaluation of the jacobian matrix)
	double kd[3][NBDIRSMAX];	// variation of the force vector at the bearing with respect to the eccentricity of the shaft e and the relative velocity of the shaft

	// kd = { { diff(k[0],e[0]), diff(k[0],e[1]), diff(k[0],e_dot[0]), diff(k[0],e_dot[1]), diff(k[0],omega_proj[0]), diff(k[0],omega_proj[1]) },
	//        { diff(k[1],e[0]), diff(k[1],e[1]), diff(k[1],e_dot[0]), diff(k[1],e_dot[1]), diff(k[1],omega_proj[0]), diff(k[1],omega_proj[1]) },
	//        { diff(k[2],e[0]), diff(k[2],e[1]), diff(k[2],e_dot[0]), diff(k[2],e_dot[1]), diff(k[2],omega_proj[0]), diff(k[2],omega_proj[1]) } };

	double eps, eps_dot, delta, SoD, SoV, my, beta;

	hydrodynamic_plain_bearing_force_dv_(m_b, m_d, m_Psi, m_eta, omega_proj, omega_projd, e_, ed, e_dot_, e_dotd, k,kd, eps,eps_dot,delta,SoD,SoV,my,beta,NBDIRSMAX);

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

	if ( m_fcontact )
	{
		// SUBROUTINE PLAIN_BEARING_CONTACT_FORCE( d, Psi, sP, DL, m, abs_FPrs1, myP, signum_delta_omega, Phi_dot, e, e_dot, k)

		plain_bearing_contact_force_dv_( m_d, m_Psi, m_sP, m_DL, m_m, m_abs_FPrs1, m_myP, m_signum_delta_omega, omega_proj, omega_projd,e_, ed, e_dot_,e_dotd, k, kd, NBDIRSMAX);

		F2_R2(1) += k[0];
		F2_R2(2) += k[1];

		dF2_R2_de_R2(1,1) += kd[0][0]; dF2_R2_de_R2(1,2) += kd[0][1];
		dF2_R2_de_R2(2,1) += kd[1][0]; dF2_R2_de_R2(2,2) += kd[1][1];

		dF2_R2_de_dot_R2(1,1) += kd[0][2]; dF2_R2_de_dot_R2(1,2) += kd[0][3];
		dF2_R2_de_dot_R2(2,1) += kd[1][2]; dF2_R2_de_dot_R2(2,2) += kd[1][3];

		dF2_R2_domega1_proj(1) += kd[0][4];
		dF2_R2_domega1_proj(2) += kd[1][4];

		dF2_R2_domega2_proj(1) += kd[0][5];
		dF2_R2_domega2_proj(2) += kd[1][5];

		M2_R2(3) += k[2];

		dM2_R2_de_R2(3,1) += kd[2][0]; 	dM2_R2_de_R2(3,2) += kd[2][1];

		dM2_R2_de_dot_R2(3,1) += kd[2][2]; 	dM2_R2_de_dot_R2(3,2) += kd[2][3];

		dM2_R2_domega1_proj(3) += kd[2][4];

		dM2_R2_domega2_proj(3) += kd[2][5];
	}

	const doublereal alpha = m_InitialAssemblyFactor.dGet();

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
	const Vec3 dlambda_dot_dX1 = -dv_dot_R2_dX1.GetRow(3) / d1_R2(3) + dv_R2_dX1.GetRow(3) / pow(d1_R2(3),2) * d1_dot_R2(3); // diff(diff(lambda,t),X1)

	const Vec3 dlambda_dX1 = -dv_R2_dX1.GetRow(3) / d1_R2(3);

	const Mat3x3 de_dot_R2_dX1 = R2.MulTM( R1.GetCol(3).Tens(dlambda_dot_dX1) + omega1.Cross( R1.GetCol(3) ).Tens(dlambda_dX1) );

	const Mat3x3 dv_dot_R2_domega1 = -R2.MulTM(Mat3x3(MatCross, R1 * m_o1_R1));
	const Mat3x3 domega1_dg1 = -Mat3x3(MatCross, omega1_ref);
	const Mat3x3 dv_R2_dg1 = -R2.MulTM( Mat3x3(MatCross, R1_0 * m_o1_R1) );
	const Mat3x3 dv_dot_R2_dg1 = dv_dot_R2_domega1 * domega1_dg1 - R2.MulTM( ( omega1 - omega2 ).Cross(Mat3x3(MatCross, R1_0 * m_o1_R1)) );
	const Mat3x3 dd1_R2_dg1 = -R2.MulTM(Mat3x3(MatCross, R1_0.GetCol(3)));
	const Mat3x3 dd1_dot_R2_dg1 = R2.MulTM( R1.GetCol(3).Cross(Mat3x3(MatCross, omega1_ref)) - ( omega1 - omega2 ).Cross(Mat3x3(MatCross, R1_0.GetCol(3))) );
	const Vec3 dlambda_dot_dg1 = -dv_dot_R2_dg1.GetRow(3) / d1_R2(3)
				    + dd1_R2_dg1.GetRow(3) * (v_dot_R2(3) / pow(d1_R2(3), 2))
				    + dv_R2_dg1.GetRow(3) / pow(d1_R2(3), 2) * d1_dot_R2(3)
				    - dd1_R2_dg1.GetRow(3) * (2. * v_R2(3) / pow(d1_R2(3), 3) * d1_dot_R2(3))
				    + dd1_dot_R2_dg1.GetRow(3) * (v_R2(3) / pow(d1_R2(3), 2));
	const Vec3 dlambda_dg1 = -dv_R2_dg1.GetRow(3) / d1_R2(3) + dd1_R2_dg1.GetRow(3) * (v_R2(3) / pow( d1_R2(3), 2 ));

	const Mat3x3 de_dot_R2_dg1 = R2.MulTM( -omega1.Cross(Mat3x3( MatCross, R1_0 * m_o1_R1 )) + ( R1 * m_o1_R1 ).Cross(Mat3x3(MatCross, omega1_ref))
								+ R1.GetCol(3).Tens(dlambda_dot_dg1) - Mat3x3( MatCross, R1_0.GetCol(3) ) * lambda_dot
								+ omega1.Cross( R1.GetCol(3) ).Tens(dlambda_dg1) + R1.GetCol(3).Cross(Mat3x3(MatCross, omega1_ref)) * lambda
								- omega1.Cross(Mat3x3( MatCross, R1_0.GetCol(3) )) * lambda );
	const Mat3x3 dv_dot_R2_dX2 = R2.MulTM(Mat3x3(MatCross, omega2));
	const Mat3x3 dv_R2_dX2 = -R2_T;
	const Vec3 dlambda_dot_dX2 = -dv_dot_R2_dX2.GetRow(3) / d1_R2(3) + dv_R2_dX2.GetRow(3) / pow(d1_R2(3),2) * d1_dot_R2(3);

	const Vec3 dlambda_dX2 = -dv_R2_dX2.GetRow(3) / d1_R2(3);
	const Mat3x3 de_dot_R2_dX2 = R2.MulTM( R1.GetCol(3).Tens(dlambda_dot_dX2) + omega1.Cross( R1.GetCol(3) ).Tens(dlambda_dX2) );
	const Mat3x3 dv_R2_dg2 = R2_0.MulTM( Mat3x3( MatCross, X1 - X2 + R1 * m_o1_R1 ) );
	const Mat3x3 dd1_R2_dg2 = R2_0.MulTM( Mat3x3( MatCross, R1 * m_o1_R1 ) );
	const Mat3x3 dv_dot_R2_dg2 = R2_0.MulTM( Mat3x3( MatCross, X1_dot - X2_dot + omega2.Cross( X2 - X1 ) + ( omega1 - omega2 ).Cross( R1 * m_o1_R1 ) ) )
				   + R2.MulTM( ( X2 - X1 - R1 * m_o1_R1 ).Cross(Mat3x3(MatCross, omega2_ref)) );
	const Vec3 dlambda_dg2 = -dv_R2_dg2.GetRow(3) / d1_R2(3) + dd1_R2_dg2.GetRow(3) * (v_R2(3) / pow(d1_R2(3),2));

	const Mat3x3 dd1_dot_R2_dg2 = R2_0.MulTM( Mat3x3( MatCross, ( omega1 - omega2 ).Cross( R1.GetCol(3) ) ) )
								- R2.MulTM( Mat3x3( MatCross, R1.GetCol(3) ) * Mat3x3(MatCross, omega2_ref) );
	const Vec3 dlambda_dot_dg2 = -dv_dot_R2_dg2.GetRow(3) / d1_R2(3) + dd1_R2_dg2.GetRow(3) * (v_dot_R2(3) / pow(d1_R2(3), 2))
				    + dv_R2_dg2.GetRow(3) / pow(d1_R2(3),2) * d1_dot_R2(3)
				    - dd1_R2_dg2.GetRow(3) * (2. * v_R2(3) / pow( d1_R2(3), 3) * d1_dot_R2(3))
				    + dd1_dot_R2_dg2.GetRow(3) * (v_R2(3) / pow(d1_R2(3),2));

	const Mat3x3 de_dot_R2_dg2 = R2_0.MulTM( Mat3x3( MatCross, X1_dot - X2_dot + omega1.Cross(R1 * m_o1_R1) - omega2.Cross( R2 * m_o2_R2 )
								+ R1.GetCol(3) * lambda_dot + omega1.Cross(R1.GetCol(3)) * lambda))
								+ R2.MulTM( -( R2 * m_o2_R2 ).Cross(Mat3x3(MatCross, omega2_ref)) + omega2.Cross(Mat3x3( MatCross, R2_0 * m_o2_R2 ))
								+ R1.GetCol(3).Tens(dlambda_dot_dg2) + omega1.Cross(R1.GetCol(3)).Tens(dlambda_dg2));
	const Mat3x3& dv_dot_R2_dX1_dot = R2_T;
	const Vec3 dlambda_dot_dX1_dot = -dv_dot_R2_dX1_dot.GetRow(3) / d1_R2(3);

	const Mat3x3 de_dot_R2_dX1_dot = R2.MulTM( Eye3 + R1.GetCol(3).Tens(dlambda_dot_dX1_dot));

	const Mat3x3 dv_dot_R2_dg1_dot = -R2.MulTM( Mat3x3( MatCross, R1 * m_o1_R1 ) );
	const Mat3x3 dd1_dot_R2_dg1_dot = -R2.MulTM( Mat3x3( MatCross, R1.GetCol(3) ) );
	const Vec3 dlambda_dot_dg1_dot = -dv_dot_R2_dg1_dot.GetRow(3) / d1_R2(3) + dd1_dot_R2_dg1_dot.GetRow(3) * (v_R2(3) / pow(d1_R2(3),2));

	const Mat3x3 de_dot_R2_dg1_dot = R2.MulTM( -Mat3x3( MatCross, R1 * m_o1_R1 ) + R1.GetCol(3).Tens(dlambda_dot_dg1_dot) - Mat3x3(MatCross, R1.GetCol(3) * lambda) );

	const Mat3x3 dv_dot_dX2_dot = -R2_T;
	const Vec3 dlambda_dot_dX2_dot = -dv_dot_dX2_dot.GetRow(3) / d1_R2(3);

	const Mat3x3 de_dot_R2_dX2_dot = R2.MulTM( -Eye3 + R1.GetCol(3).Tens(dlambda_dot_dX2_dot));

	const Mat3x3 dv_dot_R2_dg2_dot = R2.MulTM( Mat3x3( MatCross, R1 * m_o1_R1 + X1 - X2 ));
	const Mat3x3 dd1_dot_R2_dg2_dot = R2.MulTM(Mat3x3(MatCross, R1.GetCol(3)));
	const Vec3 dlambda_dot_dg2_dot = -dv_dot_R2_dg2_dot.GetRow(3) / d1_R2(3) + dd1_dot_R2_dg2_dot.GetRow(3) * (v_R2(3) / pow(d1_R2(3),2));

	const Mat3x3 de_dot_R2_dg2_dot = R2.MulTM( Mat3x3( MatCross, R2 * m_o2_R2 ) + R1.GetCol(3).Tens( dlambda_dot_dg2_dot) );

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
	WorkMat.Put(  1,  1, -dF1_I_dX1_dot - dF1_I_dX1 * dCoef );
	WorkMat.Put(  1,  4, -dF1_I_dg1_dot - dF1_I_dg1 * dCoef );
	WorkMat.Put(  1,  7, -dF1_I_dX2_dot - dF1_I_dX2 * dCoef );
	WorkMat.Put(  1, 10, -dF1_I_dg2_dot - dF1_I_dg2 * dCoef );

	WorkMat.Put(  4,  1, -dM1_I_dX1_dot - dM1_I_dX1 * dCoef );
	WorkMat.Put(  4,  4, -dM1_I_dg1_dot - dM1_I_dg1 * dCoef );
	WorkMat.Put(  4,  7, -dM1_I_dX2_dot - dM1_I_dX2 * dCoef );
	WorkMat.Put(  4, 10, -dM1_I_dg2_dot - dM1_I_dg2 * dCoef );

	WorkMat.Put(  7,  1, -dF2_I_dX1_dot - dF2_I_dX1 * dCoef );
	WorkMat.Put(  7,  4, -dF2_I_dg1_dot - dF2_I_dg1 * dCoef );
	WorkMat.Put(  7,  7, -dF2_I_dX2_dot - dF2_I_dX2 * dCoef );
	WorkMat.Put(  7, 10, -dF2_I_dg2_dot - dF2_I_dg2 * dCoef );

	WorkMat.Put( 10,  1, -dM2_I_dX1_dot - dM2_I_dX1 * dCoef );
	WorkMat.Put( 10,  4, -dM2_I_dg1_dot - dM2_I_dg1 * dCoef );
	WorkMat.Put( 10,  7, -dM2_I_dX2_dot - dM2_I_dX2 * dCoef );
	WorkMat.Put( 10, 10, -dM2_I_dg2_dot - dM2_I_dg2 * dCoef );

#ifdef DEBUG
	std::cerr << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << ":" <<  "Jac=" << std::endl << WorkMat << std::endl;
#endif
	return WorkMatV;
}




SubVectorHandler&
hydrodynamic_plain_bearing_with_offset::AssRes(SubVectorHandler& WorkVec,
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

	double eps, eps_dot, delta, SoD, SoV, my, beta;

	Vec3 e_R2, e_dot_R2;
	double omega_proj[2];
	Vec3 F2_R2, M2_R2;
	Vec3 F2_I, M2_I, F1_I, M1_I;

	ComputeResidual(e_R2,e_dot_R2,omega_proj,F2_R2,M2_R2,F2_I,M2_I,F1_I,M1_I,eps,eps_dot,delta,SoD,SoV,my,beta);
	// 1     2     3     4     5     6     7     8     9     10    11    12
	// F1(1) F1(2) F1(3) M1(1) M1(2) M1(3) F2(1) F2(2) F2(3) M2(1) M2(2) M2(3)
	WorkVec.Put(1,  F1_I);
	WorkVec.Put(4,  M1_I);
	WorkVec.Put(7,  F2_I);
	WorkVec.Put(10, M2_I);

#ifdef DEBUG
	std::cerr << __FILE__ << ":" << __LINE__ << ":" << __PRETTY_FUNCTION__ << ": Res=" << std::endl;
	std::cerr << WorkVec << std::endl;
#endif
	return WorkVec;
}

void hydrodynamic_plain_bearing_with_offset::ComputeResidual(Vec3& e_R2, Vec3& e_dot_R2,double omega_proj[2],Vec3& F2_R2,Vec3& M2_R2,Vec3& F2_I,Vec3& M2_I,Vec3& F1_I,Vec3& M1_I,double& eps, double& eps_dot,double& delta,double& SoD,double& SoV,double& my,double& beta)const
{
        const Vec3& X1 = m_pShaft->GetXCurr();
        const Vec3& X2 = m_pBearing->GetXCurr();
        const Vec3& X1_dot = m_pShaft->GetVCurr();
        const Vec3& X2_dot = m_pBearing->GetVCurr();
        const Mat3x3& R1 = m_pShaft->GetRCurr();
        const Mat3x3& R2 = m_pBearing->GetRCurr();
        const Vec3& omega1 = m_pShaft->GetWCurr();
        const Vec3& omega2 = m_pBearing->GetWCurr();

        const Vec3 v_R2 = R2.MulTV( X1 - X2 + R1 * m_o1_R1 ) - m_o2_R2;
        const Vec3 d1_R2 = R2.MulTV(R1.GetCol(3));
        const doublereal lambda = -v_R2(3) / d1_R2(3);
        e_R2 = v_R2 + d1_R2 * lambda;
        const Vec3 v_dot_R2 = R2.MulTV( X1_dot - X2_dot + omega2.Cross(X2 - X1) + ( omega1 - omega2 ).Cross( R1 * m_o1_R1 ) );
        const Vec3 d1_dot_R2 = R2.MulTV( ( omega1 - omega2 ).Cross( R1.GetCol(3) ) );
        const doublereal lambda_dot = -v_dot_R2(3) / d1_R2(3) + v_R2(3) / pow(d1_R2(3),2) * d1_dot_R2(3);

        // e_dot_R2 = R2^T * e_dot_I
        e_dot_R2 = R2.MulTV( X1_dot - X2_dot + omega1.Cross( R1 * m_o1_R1 ) - omega2.Cross( R2 * m_o2_R2 )
        					+ R1.GetCol(3) * lambda_dot + omega1.Cross(R1.GetCol(3)) * lambda);
        const Vec3 l2_R2 = m_o2_R2 + e_R2;
        const Vec3 lambda_d1_R1 = Vec3(0.,0.,lambda);
        const Vec3 l1_I = R1 * ( m_o1_R1 + lambda_d1_R1 );

        omega_proj[0] = R2.GetCol(3).Dot(omega1);
        omega_proj[1] = R2.GetCol(3).Dot(omega2);

        const double e_[2] = { e_R2(1), e_R2(2) }, e_dot_[2] = { e_dot_R2(1), e_dot_R2(2) };

        double k[3];

        hydrodynamic_plain_bearing_force_(m_b, m_d, m_Psi, m_eta, omega_proj, e_, e_dot_,k,eps,eps_dot,delta,SoD,SoV,my,beta);

        F2_R2(1) = k[0];
        F2_R2(2) = k[1];
        F2_R2(3) = 0.;

        M2_R2(1) = 0.;
        M2_R2(2) = 0.;
        M2_R2(3) = k[2];

        if ( m_fcontact )
        {
                // SUBROUTINE PLAIN_BEARING_CONTACT_FORCE( d, Psi, sP, DL, m, abs_FPrs1, myP, signum_delta_omega, Phi_dot, e, e_dot, k)

                plain_bearing_contact_force_( m_d, m_Psi, m_sP, m_DL, m_m, m_abs_FPrs1, m_myP, m_signum_delta_omega, omega_proj, e_, e_dot_, k);
                F2_R2(1) += k[0];
                F2_R2(2) += k[1];
                M2_R2(3) += k[2];
        }

        const doublereal alpha = m_InitialAssemblyFactor.dGet();

        F2_R2 *= alpha;
        M2_R2 *= alpha;

        F2_I = R2 * F2_R2;
        M2_I = R2 * ( l2_R2.Cross( F2_R2 ) + M2_R2 );
        F1_I = -F2_I;
        M1_I = -l1_I.Cross( F2_I ) - R2 * M2_R2;
}


unsigned int
hydrodynamic_plain_bearing_with_offset::iGetNumPrivData(void) const
{
	return 0;
}

int
hydrodynamic_plain_bearing_with_offset::iGetNumConnectedNodes(void) const
{
	return 2; // 1x shaft + 1x bearing
}

void
hydrodynamic_plain_bearing_with_offset::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(iGetNumConnectedNodes());
	connectedNodes[0] = m_pShaft;
	connectedNodes[1] = m_pBearing;
}

void
hydrodynamic_plain_bearing_with_offset::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{

}

std::ostream&
hydrodynamic_plain_bearing_with_offset::Restart(std::ostream& out) const
{
	out << "hydrodynamic_plain_bearing_with_offset,\n"
			"		shaft," << m_pShaft->GetLabel() << ",\n"
			"		offset," << m_o1_R1 << ",\n"
			"		bearing," << m_pBearing->GetLabel() << ",\n"
			"		offset," << m_o2_R2 << "\n"
			"		bearing_width," << m_b << ",\n"
			"		bearing_diameter," << m_d << ",\n"
			"		relative_clearance," << m_Psi << ",\n"
			"		oil_viscosity," << m_eta << ",\n"
			"		initial_assembly_factor," << m_InitialAssemblyFactor.pGetDriveCaller()->Restart(out);

	if ( m_fcontact )
	{
		out << ",\n"
			"           contact,\n"
			"			stiffness," << m_sP << ",\n"
			"			damping_ratio," << m_DL << ",\n"
			"			rated_mass," << m_m << ",\n"
			"			force," << m_abs_FPrs1 << ",\n"
			"			coulomb_friction_coefficient," << m_myP << ",\n"
			"			angular_velocity_threshold," << m_signum_delta_omega;
	}

	out << ";\n";

	return out;
}

unsigned int
hydrodynamic_plain_bearing_with_offset::iGetInitialNumDof(void) const
{
	return 0;
}

void
hydrodynamic_plain_bearing_with_offset::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
hydrodynamic_plain_bearing_with_offset::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	assert(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
hydrodynamic_plain_bearing_with_offset::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	assert(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

bool hydrodynamic_plain_bearing_set(void)
{
	UserDefinedElemRead *rf = new UDERead<hydrodynamic_plain_bearing_with_offset>;

	if (!SetUDE("hydrodynamic_plain_bearing_with_offset", rf))
	{
		delete rf;
		return false;
	}

	return true;
}

#ifndef STATIC_MODULES

#ifdef __CYGWIN__
//FIXME: dynamic linking does not work on cygwin
namespace
{	
#else
extern "C"
{
#endif

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

#ifdef __CYGWIN__
	//FIXME: dynamic linking does not work on cygwin
	volatile int g_module_init = module_init("hydrodynamic_plain_bearing",0,0);
}
#else
}
#endif

#endif // ! STATIC_MODULE

