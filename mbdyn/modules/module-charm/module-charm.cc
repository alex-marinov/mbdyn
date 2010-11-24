/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2010
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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
 * Copyright (C) 2010
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 *
 * All rights reserved.
 *
 * This module implements a user-defined induced velocity element
 * based on CDI's CHARM free wake.
 *
 * <http://www.continuum-dynamics.com/pr-charm.html>
 *
 * This module has been sponsored by Baldwin Technology Company LLC
 * <http://www.baldwintechnology.com/>
 */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"
#include "indvel.h"
#include "drive_.h"

// CHARM's header file (apparently it is experimental)
#undef __stdcall
#define __stdcall
#include <CharmWP.h>

// The functions that follow are stubs to simulate the availability
// of CHARM's WP library
#define MBCHARM_FAKE 1
#ifdef MBCHARM_FAKE
int
allocWPAircraft(
	int num_rotor_surf,
	int *num_blades,
	int *num_control_pts,
	wpAircraft *aircraft)
{
	ASSERT(num_rotor_surf > 0);
	ASSERT(num_blades != 0);
	ASSERT(num_control_pts != 0);
	ASSERT(aircraft != 0);

	aircraft->number_rotors_or_surfaces = num_rotor_surf;

	for (int ir = 0; ir < num_rotor_surf; ir++) {
		wpRotorSurface *rotor = &aircraft->rotors[ir];

		ASSERT(num_blades[ir] > 0);
		ASSERT(num_control_pts[ir] > 0);

		rotor->number_blades = num_blades[ir];
		rotor->number_control_points = num_control_pts[ir];
	}

	return 0;
}

void
initWPModule( 
	int num_aircraft, 
	wpAircraft *aircraft, 
	int num_offrotor, 
	int max_offrotor, 
	chOffRotorEval *eval, 
	int areEvalPointsFixed, 
	chGlobal *global, 
	wpOptions *options)
{
	ASSERT(num_aircraft == 1);
	ASSERT(aircraft != 0);
	ASSERT((num_offrotor == 0) || (eval != 0));
	ASSERT(global != 0);
	ASSERT(options != 0);
}

void
updateIndVelocity(wpAircraft *aircraft, chOffRotorEval *eval)
{
	ASSERT(aircraft != 0);
	ASSERT(eval != 0);
}

void
setCHARMGlobal(chGlobal *global)
{
	ASSERT(global != 0);
}

void
setWPAircraft(int index, wpAircraft *aircraft)
{
	ASSERT(aircraft != NULL);
	ASSERT(index == 0);	// only one aircraft
}

void
setCHARMOffRotor(int numPoints, chOffRotorEval *eval)
{
	ASSERT(eval != 0);
}

void
updateWake(int isTrimming, chStatus *status)
{
	ASSERT(status != 0);
}
#endif // MBCHARM_FAKE

class ModuleCHARM
: virtual public Elem,
	public UserDefinedElem,
	public InducedVelocity
{
private:
	// TODO: define per-point structure
	struct PointData {
		unsigned label;
		unsigned counter;
		Vec3 F;
		Vec3 M;
		doublereal dW;
		Vec3 X;

	public:
		PointData(void) { NO_OP; };
		~PointData(void) { NO_OP; };
	};

	std::vector<PointData> m_data;
	typedef std::vector<PointData> PD;
	PD::iterator m_data_iter;

	// rotor->blade->point
	struct ElemMapping {
		unsigned uLabel;
	};

	struct SurfaceMapping {
		std::vector<ElemMapping> Elems;
	};

	struct RotorMapping {
		std::vector<SurfaceMapping> Blades;
	};

	std::vector<RotorMapping> m_Rotors;

	// element -> rotor, blade
	struct RotorBlade {
		int iRotor;
		int iBlade;
		int iElem;
		int iFirstPoint;

		RotorBlade(int iRotor, int iBlade, int iElem)
			: iRotor(iRotor), iBlade(iBlade),
			iElem(iElem), iFirstPoint(-1)
		{
			NO_OP;
		};

		virtual ~RotorBlade(void) { NO_OP; };
	};

	typedef std::map<unsigned, RotorBlade *> RBM;
	RBM m_e2b_map;

	// Time handling
	// KTRSIM-related stuff (Fidelity/speed trade off settings)
	// TODO: use azimuth instead of time (number of revolutions)
	DriveOwner m_Time;
	std::vector<double> m_TrimTime;
	std::vector<double>::const_iterator m_TrimTimeIter;

	chGlobal m_chglobal;
	wpOptions m_wpoptions;
	std::vector<chOffRotorEval> m_eval_pts;

	Mat3x3 m_Rh;
	wpAircraft m_wpaircraft;

	// add private data
	int iFirstAssembly;

	void Init_int(void);
	void Update_int(void);

public:
	ModuleCHARM(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleCHARM(void);

	virtual Elem::Type GetElemType(void) const;

	// induced velocity specific calls
	virtual InducedVelocity::Type GetInducedVelocityType(void) const;
	virtual bool bSectionalForces(void) const;
	virtual Vec3 GetInducedVelocity(const Vec3&) const;
	virtual void AddSectionalForce(unsigned int uL, unsigned iPnt,
		const Vec3& F, const Vec3& M, doublereal dW,
		const Vec3& X, const Mat3x3& R,
		const Vec3& V, const Vec3& W);

	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
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
};

ModuleCHARM::ModuleCHARM(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
InducedVelocity(uLabel, 0, 0, flag(0)),
m_Time(new TimeDriveCaller(pDM->pGetDrvHdl())),
iFirstAssembly(2)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	CHARM							\n"
"Author: 	Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
" Description:	This module implements an induced velocity model	\n"
"		based on CHARM's free wake, distributed by CDI.		\n"
"		Sponsored by Baldwin Technologies Company, LLC.		\n"
"									\n"
"	All rights reserved.						\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// global data
	m_chglobal.current_time = m_Time.dGet();
	m_chglobal.trim_flag = 0;

	// NOTE: ignore ground by now
	m_chglobal.ground_flag = 0;

	// NOTE: ignore filters
	m_chglobal.iv_filter_order[0] = 0;
	m_chglobal.iv_filter_order[1] = 0;
#if 0
	m_chglobal.iv_filter_freq[0];
	m_chglobal.iv_filter_freq[1];
	m_chglobal.iv_filter_damp[0];
	m_chglobal.iv_filter_damp[1];
#endif

	// wake-panel options
	m_wpoptions.kblind = 1;
	m_wpoptions.klssim = 0;
	m_wpoptions.kgmsim = 0;
	m_wpoptions.kupdgam = 1;	// set to 0 when KTRSIM is set to 4
	m_wpoptions.kshipsim = 0;
	m_wpoptions.kc2sim = 0;
	m_wpoptions.kc2lev = 0;
	if (!HP.IsKeyWord("units")) {
		silent_cerr("ModuleCHARM(" << uLabel << "): "
			"\"units\" expected "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("bu")) {
		m_wpoptions.kunitsim = 0;	// 0: BU; 1: MKS

	} else if (HP.IsKeyWord("si") || HP.IsKeyWord("mks")) {
		m_wpoptions.kunitsim = 1;	// 0: BU; 1: MKS

	} else {
		silent_cerr("ModuleCHARM(" << uLabel << "): "
			"unknown \"units\" "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_wpoptions.noutsim = 998;
	if (HP.IsKeyWord("output" "unit")) {
		m_wpoptions.noutsim = HP.GetInt();
		if (m_wpoptions.noutsim <= 0) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"invalid output unit " << m_wpoptions.noutsim << " "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	m_wpoptions.nupsim = 0;


	// aircraft
	if (!HP.IsKeyWord("aircraft")) {
		silent_cerr("ModuleCHARM(" << uLabel << "): "
			"\"aircraft\" (aircraft node label) expected "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pCraft = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

	ReferenceFrame RF(pCraft);
	if (HP.IsKeyWord("orientation")) {
		m_Rh = HP.GetRotRel(RF);
	} else {
		m_Rh = Eye3;
	}

	if (!HP.IsKeyWord("rotors")) {
		silent_cerr("ModuleCHARM(" << uLabel << "): "
			"\"rotors\" (number of rotors) expected "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	int num_rotors = HP.GetInt();
	if (num_rotors < 1) {
		silent_cerr("ModuleCHARM(" << uLabel << "): "
			"invalid number of rotors " << num_rotors
			<< " at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_Rotors.resize(num_rotors);

	for (int ir = 0; ir < num_rotors; ir++) {
		int iNumBlades = HP.GetInt();
		if (iNumBlades < 1) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"invalid number of blades " << iNumBlades
				<< " for rotor #" << ir << " of " << num_rotors
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_Rotors[ir].Blades.resize(iNumBlades);

		int iNumElems = HP.GetInt();
		if (iNumElems < 1) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"invalid number of elements " << iNumElems
				<< " for rotor #" << ir << " of " << num_rotors
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		for (int ib = 0; ib < iNumBlades; ib++) {
			m_Rotors[ir].Blades[ib].Elems.resize(iNumElems);
			for (int ie = 0; ie < iNumElems; ie++) {
				int il = HP.GetInt();
				if (il < 0) {
					silent_cerr("ModuleCHARM(" << uLabel << "): "
						"invalid label " << il
						<< " for element #" << ie << " of " << iNumElems << ", "
						<< "blade #" << ib << " of " << iNumBlades << ", "
						<< "rotor #" << ir << " of " << num_rotors
						<< " at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				RotorBlade *pRB = new RotorBlade(ir, ib, ie);
				if (!m_e2b_map.insert(RBM::value_type(unsigned(il), pRB)).second) {
					delete pRB;
					silent_cerr("ModuleCHARM(" << uLabel << "): "
						"label " << il
						<< " for element #" << ie << " of " << iNumElems << ", "
						<< "blade #" << ib << " of " << iNumBlades << ", "
						<< "rotor #" << ir << " of " << num_rotors
						<< "already in use"
						<< " at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				m_Rotors[ir].Blades[ib].Elems[ie].uLabel = unsigned(il);
			}
		}
	}

	if (HP.IsKeyWord("trim" "times")) {
		m_TrimTime.resize(5);
		for (unsigned it = 0; it < 4; it++) {
			m_TrimTime[it] = HP.GetReal();
			if (it > 0 && m_TrimTime[it] <= m_TrimTime[it - 1]) {
				silent_cerr("ModuleCHARM(" << uLabel << "): "
					"TrimTime[" << it << "]=" << m_TrimTime[it] << " "
					"incompatible with TrimTime[" << it - 1 << "]=" << m_TrimTime[it - 1] << " "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

	} else {
		m_TrimTime.resize(1);
	}
	m_TrimTime[m_TrimTime.size() - 1] = std::numeric_limits<doublereal>::max();
	m_TrimTimeIter = m_TrimTime.begin();

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

ModuleCHARM::~ModuleCHARM(void)
{
	// destroy private data
	NO_OP;
}

Elem::Type
ModuleCHARM::GetElemType(void) const
{
	return Elem::LOADABLE;
}

InducedVelocity::Type
ModuleCHARM::GetInducedVelocityType(void) const
{
	return InducedVelocity::USER_DEFINED;
}

bool
ModuleCHARM::bSectionalForces(void) const
{
	return true;
}

void
ModuleCHARM::Init_int(void)
{
	// this function must be called after all data has been set up
	ASSERT(iFirstAssembly == 1);

	int num_rotors = m_Rotors.size();
	std::vector<int> v_num_blades(num_rotors);
	std::vector<int> v_num_control_points(num_rotors);

	for (int ir = 0; ir < num_rotors; ir++) {
		int num_blades = m_Rotors[ir].Blades.size();
		ASSERT(num_blades > 0);
		int num_points_per_blade = 0;
		for (int ib = 0; ib < num_blades; ib++) {
			int num_elems = m_Rotors[ir].Blades[ib].Elems.size();
			int num_points = 0;
			for (int ie = 0; ie < num_elems; ie++) {
				unsigned uLabel = m_Rotors[ir].Blades[ib].Elems[ie].uLabel;
				int iPoints = m_e2b_map[uLabel]->iFirstPoint;

#if 0
				std::cerr << "ModuleCHARM(" << GetLabel() << ")::Init_int: "
					"label=" << uLabel
					<< " ir=" << m_e2b_map[uLabel]->iRotor
					<< " ib=" << m_e2b_map[uLabel]->iBlade
					<< " ie=" << m_e2b_map[uLabel]->iElem
					<< " num_points=" << m_e2b_map[uLabel]->iFirstPoint
					<< std::endl;
#endif

				if (iPoints == -1) {
					silent_cerr("ModuleCHARM(" << uLabel << "): "
						"Element(" << uLabel << "), "
						"element " << ie << " of " << num_elems << ", "
						"blade " << ib << " of " << num_blades << ", "
						"rotor " << ir << " of " << num_rotors << " "
						"did not write forces during initialization" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				num_points += iPoints;
			}

#if 0
			std::cerr << "ModuleCHARM(" << uLabel << "): "
				"rotor=" << ir << " "
				"blade=" << ib << " "
				"points=" << num_points << std::endl;
#endif

			if (ib == 0) {
				num_points_per_blade = num_points;

			} else if (num_points != num_points_per_blade) {
				silent_cerr("ModuleCHARM(" << uLabel << "): "
					"number of points " << num_points << ", "
					"blade " << ib << " of " << num_blades << ", "
					"rotor " << ir << " of " << num_rotors << " "
					"does not match number of points " << num_points_per_blade << " "
					"of blade " << 0 << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		v_num_control_points[ir] = num_points_per_blade;
		if (v_num_control_points[ir] < 1) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"invalid number of control points " << v_num_control_points[ir]
				<< " for rotor #" << ir << " of " << num_rotors << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	allocWPAircraft(num_rotors, &v_num_blades[0],
           	&v_num_control_points[0], &m_wpaircraft);

	if (m_wpaircraft.number_rotors_or_surfaces != num_rotors) {
		silent_cerr("ModuleCHARM(" << uLabel << "): "
			"number of rotors mismatch "
			"(" << m_wpaircraft.number_rotors_or_surfaces
			<< " instead of " << num_rotors << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (int ir = 0; ir < num_rotors; ir++) {
		wpRotorSurface *rotor = &m_wpaircraft.rotors[ir];
		if (rotor->number_blades != v_num_blades[ir]) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"number of blades mismatch "
				"(" << rotor->number_blades
				<< " instead of " << v_num_blades[ir] << ") "
				"for rotor #" << ir << " of " << num_rotors << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (rotor->number_control_points != v_num_control_points[ir]) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"number of control points mismatch "
				"(" << rotor->number_control_points
				<< " instead of " << v_num_control_points[ir] << ") "
				"for rotor #" << ir << " of " << num_rotors << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	// FIXME: one aircraft only
	initWPModule(1,
		&m_wpaircraft,
		m_eval_pts.size(), m_eval_pts.size(), &m_eval_pts[0],
		1,		// FIXME: areEvalPointsFixed?
		&m_chglobal, &m_wpoptions);
}

void
ModuleCHARM::Update_int(void)
{
	const Vec3& Xac(pCraft->GetXCurr());
	// p, T not used yet
	doublereal rho, c, p, T;
	GetAirProps(Xac, rho, c, p, T);
	const Vec3& Vac(pCraft->GetVCurr());
	Mat3x3 Rac(pCraft->GetRCurr()*m_Rh);
	Vec3 VBac(Rac.MulTV(Vac));
	Vec3 WBac(Rac.MulTV(pCraft->GetWCurr()));
	// FIXME: Euler angles? 123 or 321?  Rac or Rac^T?
	Vec3 EBac = MatR2EulerAngles123(Rac);

	// update global data
	Vec3 Vinf(0.);
	if (fGetAirVelocity(Vinf, Xac)) {
		m_chglobal.inertial_winds[0] = Vinf(1);
		m_chglobal.inertial_winds[1] = Vinf(2);
		m_chglobal.inertial_winds[2] = Vinf(3);
	}

	doublereal dTime = m_Time.dGet();
	m_chglobal.time_increment = dTime - m_chglobal.current_time;
	m_chglobal.current_time = dTime;

	// NOTE: ignore ground by now
#if 0
	m_chglobal.ground_position[0];
	m_chglobal.ground_position[1];
	m_chglobal.ground_position[2];
	m_chglobal.ground_orientation[0];
	m_chglobal.ground_orientation[1];
	m_chglobal.ground_orientation[2];
	m_chglobal.ground_center[0];
	m_chglobal.ground_center[1];
	m_chglobal.ground_center[2];
#endif

	// time-driven trim flag
	while (dTime >= *m_TrimTimeIter) {
		m_TrimTimeIter++;
		m_chglobal.trim_flag++;
		if (m_chglobal.trim_flag == 4) {
			// FIXME: kupdgam may need to remain 1
			m_wpoptions.kupdgam = 0;
		}
		silent_cout("ModuleCHARM(" << GetLabel() << "): "
			"Time=" << dTime << " trim_flag=" << m_chglobal.trim_flag << std::endl);
	}

	// update aircraft data
	m_wpaircraft.air_density = rho;
	m_wpaircraft.inertial_position[0] = Xac(1);
	m_wpaircraft.inertial_position[1] = Xac(2);
	m_wpaircraft.inertial_position[2] = Xac(3);
	m_wpaircraft.inertial_velocity[0] = Vac(1);
	m_wpaircraft.inertial_velocity[1] = Vac(2);
	m_wpaircraft.inertial_velocity[2] = Vac(3);
	m_wpaircraft.body_velocity[0] = VBac(1);
	m_wpaircraft.body_velocity[1] = VBac(2);
	m_wpaircraft.body_velocity[2] = VBac(3);
	m_wpaircraft.angular_rates[0] = WBac(1);
	m_wpaircraft.angular_rates[1] = WBac(2);
	m_wpaircraft.angular_rates[2] = WBac(3);
	m_wpaircraft.euler_angles[0] = EBac(1);
	m_wpaircraft.euler_angles[1] = EBac(2);
	m_wpaircraft.euler_angles[2] = EBac(3);
	// FIXME: R(i, j) or R(j, i)?
	m_wpaircraft.T_inertial_to_body[0][0] = Rac(1, 1);
	m_wpaircraft.T_inertial_to_body[0][1] = Rac(2, 1);
	m_wpaircraft.T_inertial_to_body[0][2] = Rac(3, 1);
	m_wpaircraft.T_inertial_to_body[1][0] = Rac(1, 2);
	m_wpaircraft.T_inertial_to_body[1][1] = Rac(2, 2);
	m_wpaircraft.T_inertial_to_body[1][2] = Rac(3, 2);
	m_wpaircraft.T_inertial_to_body[2][0] = Rac(1, 3);
	m_wpaircraft.T_inertial_to_body[2][1] = Rac(2, 3);
	m_wpaircraft.T_inertial_to_body[2][2] = Rac(3, 3);

	// TODO: rotor(s)

	// actually call the module and update the induced velocity
	// FIXME: the exact sequence needs to be carefully defined

	// evaluate induced velocity at rotor blade evaluation points
	// and off-rotor evaluation points
	updateIndVelocity(&m_wpaircraft, &m_eval_pts[0]);

	// update if needed
	setCHARMGlobal(&m_chglobal);

	// update if needed
	setWPAircraft(0, &m_wpaircraft);

	// update if needed
	setCHARMOffRotor(m_eval_pts.size(), &m_eval_pts[0]);

	// update wake
	int isTrimming = 0;
	chStatus status;
	updateWake(isTrimming, &status);
	// FIXME: test status?
}

Vec3
ModuleCHARM::GetInducedVelocity(const Vec3& X) const
{
	return Zero3;
}

void
ModuleCHARM::AddSectionalForce(unsigned int uL, unsigned iPnt,
	const Vec3& F, const Vec3& M, doublereal dW,
	const Vec3& X, const Mat3x3& R, const Vec3& V, const Vec3& W)
{
	std::cerr << "ModuleCHARM(" << GetLabel() << ")::AddSectionalForce: " << uL << ":" << iPnt << std::endl;

	if (iFirstAssembly == 1) {
		unsigned idx = m_data.size();
		m_data.resize(idx + 1);
		m_data[idx].label = uL;
		if (idx == 0 || (idx > 0 && uL != m_data[idx - 1].label)) {
			m_data[idx].counter = 0;
		} else {
			m_data[idx].counter = m_data[idx - 1].counter + 1;
		}

		m_data[idx].F = Zero3;
		m_data[idx].M = Zero3;
		m_data[idx].dW = 0.;
		m_data[idx].X = Zero3;

		int ir = m_e2b_map[uL]->iRotor;
		int ib = m_e2b_map[uL]->iBlade;
		int ie = m_e2b_map[uL]->iElem;
		if (m_Rotors[ir].Blades[ib].Elems[ie].uLabel != uL) {
			silent_cerr("ModuleCHARM(" << GetLabel() << "): "
				"uL=" << uL << " iPnt=" << iPnt << " mismatch" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (m_e2b_map[uL]->iFirstPoint == -1) {
			m_e2b_map[uL]->iFirstPoint = 1;
		} else {
			m_e2b_map[uL]->iFirstPoint++;
		}

#if 0
		std::cerr << "ModuleCHARM(" << GetLabel() << ")::AddSectionalForce: " << uL << ":" << iPnt
			<< " ir=" << ir << " ib=" << ib << " ie=" << ie << " num_points=" << m_e2b_map[uL]->iFirstPoint
			<< std::endl;
#endif

	} else {
		Mat3x3 R(pCraft->GetRCurr()*m_Rh);

		// resolve force, moment and point in craft's reference frame
		m_data_iter->F = R.MulTV(F);
		m_data_iter->M = R.MulTV(M);
		m_data_iter->dW = dW;
		m_data_iter->X = R.MulTV(X - pCraft->GetXCurr());

		m_data_iter++;
	}
}

void
ModuleCHARM::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	// NO_OP;
	Update_int();
}

void
ModuleCHARM::Output(OutputHandler& OH) const
{
	std::cerr << "ModuleCHARM(" << GetLabel() << ")::Output: size=" << m_data.size() << std::endl;

	// should do something useful
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();

		for (PD::const_iterator i = m_data.begin(); i != m_data.end(); i++) {
			out << GetLabel() << "#" << i->label << "#" << i->counter
				<< " " << i->F
				<< " " << i->M
				<< " " << i->dW
				<< " " << i->X
				<< std::endl;
		}
	}
}

void
ModuleCHARM::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler& 
ModuleCHARM::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleCHARM::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	std::cerr << "ModuleCHARM(" << GetLabel() << ")::AssRes: iFirstAssembly=" << iFirstAssembly << std::endl;

	if (iFirstAssembly) {
		if (iFirstAssembly == 1) {
			// all topological information should be available;
			// module can be initialized
			Init_int();
		}

		iFirstAssembly--;
	}

	// re-initialize iterator to loop over point data
	// when AddSectionalForces is called
	m_data_iter = m_data.begin();

	// should do something useful
	WorkVec.ResizeReset(0);

	return WorkVec;
}

unsigned int
ModuleCHARM::iGetNumPrivData(void) const
{
	return 0;
}

int
ModuleCHARM::iGetNumConnectedNodes(void) const
{
	return 1;
}

void
ModuleCHARM::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);
	connectedNodes[0] = pCraft;
}

void
ModuleCHARM::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ModuleCHARM::Restart(std::ostream& out) const
{
	return out << "# ModuleCHARM: not implemented" << std::endl;
}

unsigned int
ModuleCHARM::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleCHARM::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleCHARM::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleCHARM::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<ModuleCHARM>;

	if (!SetUDE("charm", rf)) {
		delete rf;

		silent_cerr("module-charm: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

