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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstdlib>
#include <cerrno>
#include <cfloat>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "dataman.h"
#include "userelem.h"
#include "indvel.h"
#include "drive_.h" // for TimeDrive
#include "Rot.hh"

#include "mbcharm.h"

// CHARM's header file (apparently it is experimental)
#undef __stdcall
#define __stdcall
#include <CharmWP.h>

// The functions that follow are stubs to simulate the availability
// of CHARM's WP library
//#define MBCHARM_FAKE 1
#ifdef MBCHARM_FAKE
extern "C" int
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

extern "C" void
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

extern "C" void
updateIndVelocity(wpAircraft *aircraft, chOffRotorEval *eval)
{
	ASSERT(aircraft != 0);
	ASSERT(eval != 0);
}

extern "C" void
setCHARMGlobal(chGlobal *global)
{
	ASSERT(global != 0);
}

extern "C" void
setWPAircraft(int index, wpAircraft *aircraft)
{
	ASSERT(aircraft != NULL);
	ASSERT(index == 0);	// only one aircraft
}

extern "C" void
setCHARMOffRotor(int numPoints, chOffRotorEval *eval)
{
	ASSERT(eval != 0);
}

extern "C" void
updateWake(int isTrimming, chStatus *status)
{
	ASSERT(status != 0);
}
#endif // MBCHARM_FAKE

#ifdef NEED_CHARM_SHIPWAKE
// Some builds of libWPModule call this function...

/*
      SUBROUTINE SHIPWAKE(T,NPTS,XPTS,VPTS,IFRAME)
C   This subroutine provides ship airwake induced velocities
C   at a set of points at time T.  Ship airwake velocities
C   are added to the current values in the VPTS array.

C--------------------C
C  List of Arguments C
C--------------------C
C  
C   ARGUMENT       TYPE       I/O   MEANING
C
C    T             REAL        I    Time (seconds)
C
C    NPTS          INT         I    Number of evaluation points where
C                                   ship airwake induced velocities
C                                   are to be evaluated at time T
C
C    XPTS(i,n)  REAL(3,NPTS)   I    Evaluation point locations at time T
C                                   i=1,2,3 = x,y,z; n=point number
C                                    (ft or m)
C
C    VPTS(i,n)  REAL(3,NPTS)   O    The ship airwake induced velocity
C                                   at the XPTS evaluation points at
C                                   time T is added to the velocity
C                                   currently in VPTS. (ft/s or m/s)
C
C    IFRAME         INT        I    IFRAME=1 points and velocities are 
C                                            in the ship frame
C                                   IFRAME=2 points and velocities are
C                                            in the inertial frame
C--------------------C
C  List of Variables C
C--------------------C
C  
C    XSHIP(i)      REAL(3)     Origin of ship airwake frame
C                              in inertial frame at time T,
C                              (i=1,2,3 = x,y,z) (ft or m)
C  
C    VSHIP(i)      REAL(3)     Velocity of ship airwake frame
C                              in inertial frame at time T,
C                              (i=1,2,3 = u,v,w) (ft/s or m/s)
C  
C    ROLL,PITCH    REAL        Orientation of ship airwake frame
C    YAW                       in inertial frame at time T,
C                              (degs) (applied in reverse order)
C
C    WINDMAG       REAL        Incoming free stream wind magnitude 
C                               (ft/s or m/s)
C
C    WINDDIR       REAL        Incoming free stream wind direction
C                                (degs) (positive from starboard)
C
C    WIND(i)       REAL(3)     Wind vector in the ship airwake frame
C                              at time T, (i=1,2,3 = u,v,w)
C                              (ft/s or m/s) (determined from
C                              WINDMAG and WINDDIR).  WIND is the 
C                              incoming wind vector.
C
C    ITEST         INT         Flags the type of ship airwake model
C                              ITEST = 0: Use the ship airwake
C                                         read in from the shipwake.dat
C                                         file
C                              ITEST > 0: Use a user-defined analytical
C                                         ship airwake model coded up
C                                         directly in this subroutine.

C---------------------------------------------------------------------C
C** NOTE: The ship airwake induced velocity is *added* to VPTS,
C**       as opposed to simply assigned to VPTS.  This saves 
C**       memory and some CPU time.
C---------------------------------------------------------------------C

      DIMENSION XSHIP(3),VSHIP(3),WIND(3),XPTS(3,*),VPTS(3,*)
      DIMENSION TV(3),XV(3),VV(3),FITS(3,3)
      DATA ITEST/1/

      DATA RTD,DTR/57.2957795,.017453293/
 */

extern "C" int
shipwake_(
	doublereal *T,
	integer *NPTS,
	doublereal *XPTS,
	doublereal *VPTS,
	integer *IFRAME )
{
	return 0;
}

#endif // ! NEED_CHARM_SHIPWAKE

class ModuleCHARM
: virtual public Elem,
	public UserDefinedElem,
	public InducedVelocity
{
private:
	// save data manager
	DataManager *pDM;

	// forward declaration
	struct RotorBlade;

	// TODO: define per-point structure
	struct PointData {
		Elem::Type type;
		unsigned label;
		unsigned counter;

		RotorBlade *pRB;
		int iOff;

		Vec3 X;

		double *pos;
		double *vel;

		double spanwise_lift;
		double tangential_velocity;

		double *spanwise_lift_p;
		double *tangential_velocity_p;

	public:
		PointData(void)
		: type(Elem::UNKNOWN), label(unsigned(-1)), counter(unsigned(-1)),
		pRB(0), iOff(-1), pos(0), vel(0),
		spanwise_lift_p(0), tangential_velocity_p(0)
		{
			NO_OP;
		};
		~PointData(void) { NO_OP; };
	};

	mutable std::vector<PointData> m_data;
	typedef std::vector<PointData> PD;
	PD::iterator m_data_frc_iter;
	mutable PD::iterator m_data_vel_iter;

	// element -> rotor, blade
	struct RotorBlade {
		int iRotor;
		int iBlade;
		int iElem;
		int iFirst;
		int iCount;
		int iIdx;

		RotorBlade(int iRotor, int iBlade, int iElem)
			: iRotor(iRotor), iBlade(iBlade),
			iElem(iElem), iCount(-1), iIdx(-1)
		{
			NO_OP;
		};

		virtual ~RotorBlade(void) { NO_OP; };
	};

	// TODO: need to take element type into account
	typedef std::map<unsigned, RotorBlade *> RBM;
	RBM m_e2b_map;

	// rotor->blade->point
	struct ElemMapping {
		unsigned uLabel;
		RotorBlade *pRB;
	};

	struct SurfaceMapping {
		Mat3x3 Rh;

		std::vector<ElemMapping> Elems;
	};

	struct RotorMapping {
		int rotation_dir;
		double radius;
		double average_chord;
		double root_cutout;
		double omega100;

		Mat3x3 Rh;
		StructNode *pHub;

		std::vector<SurfaceMapping> Blades;
	};

	std::vector<RotorMapping> m_Rotors;

	// Time handling
	// KTRSIM-related stuff (Fidelity/speed trade off settings)
	// TODO: use azimuth instead of time (number of revolutions)
	std::vector<double> m_TrimTime;
	std::vector<double>::const_iterator m_TrimTimeIter;

	chGlobal m_chglobal;
	wpOptions m_wpoptions;
	std::vector<chOffRotorEval> m_eval_pts;

	Mat3x3 m_Rh;
	Mat3x3 m_Rac;
	wpAircraft m_wpaircraft;

	// add private data
	int iFirstAssembly;

	void Init_int(void);
	void Set_int(void);
	void Update_int(void);

public:
	ModuleCHARM(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleCHARM(void);

	virtual Elem::Type GetElemType(void) const;

	// induced velocity specific calls
	virtual InducedVelocity::Type GetInducedVelocityType(void) const;
	virtual bool bSectionalForces(void) const;
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3&) const;
	virtual void AddSectionalForce(Elem::Type type,
		unsigned int uLabel, unsigned uPnt,
		const Vec3& F, const Vec3& M, doublereal dW,
		const Vec3& X, const Mat3x3& R,
		const Vec3& V, const Vec3& W);

	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);
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
pDM(pDM),
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

	m_chglobal.inertial_winds[0] = 0.;
	m_chglobal.inertial_winds[1] = 0.;
	m_chglobal.inertial_winds[2] = 0.;

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
		if (HP.IsKeyWord("rotation" "direction")) {
			if (HP.IsKeyWord("counter" "clockwise")) {
				m_Rotors[ir].rotation_dir = 1;

			} else if (HP.IsKeyWord("clockwise")) {
				m_Rotors[ir].rotation_dir = -1;

			} else {
				m_Rotors[ir].rotation_dir = HP.GetInt();
				if (std::abs(m_Rotors[ir].rotation_dir) != 1) {
					silent_cerr("ModuleCHARM(" << uLabel << "): "
						"invalid \"rotation direction\" "
						"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}

		} else {
			m_Rotors[ir].rotation_dir = 1;
		}

		if (!HP.IsKeyWord("radius")) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"\"radius\" expected "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_Rotors[ir].radius = HP.GetReal();
		if (m_Rotors[ir].radius <= 0.) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"invalid \"radius\" " << m_Rotors[ir].radius
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("chord")) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"\"chord\" expected "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_Rotors[ir].average_chord = HP.GetReal();
		if (m_Rotors[ir].average_chord <= 0.) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"invalid (average) \"chord\" " << m_Rotors[ir].average_chord
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("root" "cutout")) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"\"root cutout\" expected "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_Rotors[ir].root_cutout = HP.GetReal();
		if (m_Rotors[ir].root_cutout <= 0.
			|| m_Rotors[ir].root_cutout >= m_Rotors[ir].radius)
		{
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"invalid \"root cutout\" " << m_Rotors[ir].root_cutout
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("omega")) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"\"omega\" expected "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_Rotors[ir].omega100 = HP.GetReal();
		if (m_Rotors[ir].omega100 <= 0.) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"invalid \"root cutout\" " << m_Rotors[ir].omega100
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("hub" "node")) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"\"hub node\" expected "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_Rotors[ir].pHub = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

		if (HP.IsKeyWord("orientation")) {
			m_Rotors[ir].Rh = HP.GetRotRel(ReferenceFrame(pCraft));

		} else {
			// construct Rh from aircraft and hub nodes
			Mat3x3 Rac(pCraft->GetRCurr()*m_Rh);
			m_Rotors[ir].Rh = Rac.MulTM(m_Rotors[ir].pHub->GetRCurr());
		}

		if (!HP.IsKeyWord("blades")) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"\"blades\" expected "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		int iNumBlades = HP.GetInt();
		if (iNumBlades < 1) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"invalid number of blades " << iNumBlades
				<< " for rotor #" << ir << " of " << num_rotors
				<< " at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		m_Rotors[ir].Blades.resize(iNumBlades);

		if (!HP.IsKeyWord("elements")) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"\"elements\" expected "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
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
				m_Rotors[ir].Blades[ib].Elems[ie].pRB = pRB;
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

#if 0
	std::cerr << "ModuleCHARM(" << GetLabel() << ")::Init_int: "
		"num_rotors=" << num_rotors << std::endl;
#endif

	for (int ir = 0; ir < num_rotors; ir++) {
		int num_blades = m_Rotors[ir].Blades.size();

#if 0
		std::cerr << "ModuleCHARM(" << GetLabel() << ")::Init_int: "
			"rotor[" << ir << "] num_blades=" << num_blades << std::endl;
#endif

		ASSERT(num_blades > 0);
		int num_points_per_blade = 0;
		for (int ib = 0; ib < num_blades; ib++) {
			int num_elems = m_Rotors[ir].Blades[ib].Elems.size();
			int num_points = 0;
			for (int ie = 0; ie < num_elems; ie++) {
				unsigned uLabel = m_Rotors[ir].Blades[ib].Elems[ie].uLabel;
				int iPoints = m_e2b_map[uLabel]->iCount;

#if 0
				std::cerr << "ModuleCHARM(" << GetLabel() << ")::Init_int: "
					"label=" << uLabel
					<< " ir=" << m_e2b_map[uLabel]->iRotor
					<< " ib=" << m_e2b_map[uLabel]->iBlade
					<< " ie=" << m_e2b_map[uLabel]->iElem
					<< " num_points=" << m_e2b_map[uLabel]->iCount
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

				m_e2b_map[uLabel]->iFirst = num_points;
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

		if (num_points_per_blade < 1) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"invalid number of control points " << v_num_control_points[ir]
				<< " for rotor #" << ir << " of " << num_rotors << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		v_num_blades[ir] = num_blades;
		v_num_control_points[ir] = num_points_per_blade;
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

#if 0
	std::cerr << "ModuleCHARM(" << GetLabel() << ")::Init_int: "
		"num_rotors_or_surfaces=" << m_wpaircraft.number_rotors_or_surfaces << std::endl;
#endif

	Mat3x3 Rb2h(1, 0, 0, 0, -1, 0, 0, 0, -1);

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

#if 0
		std::cerr << "ModuleCHARM(" << GetLabel() << ")::Init_int: "
			"rotor[" << ir << "] number_blades=" << rotor->number_blades << std::endl;
#endif

		if (rotor->number_control_points != v_num_control_points[ir]) {
			silent_cerr("ModuleCHARM(" << uLabel << "): "
				"number of control points mismatch "
				"(" << rotor->number_control_points
				<< " instead of " << v_num_control_points[ir] << ") "
				"for rotor #" << ir << " of " << num_rotors << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		rotor->isRotor = 1;
		ASSERT(std::abs(m_Rotors[ir].rotation_dir) == 1);
		rotor->rotation_dir = m_Rotors[ir].rotation_dir;
		rotor->radius = m_Rotors[ir].radius;
		rotor->average_chord = m_Rotors[ir].average_chord;
		rotor->root_cutout = m_Rotors[ir].root_cutout;
		// T_body_to_hub = Rh^T
		// NOTE: C indexes need to be exchanged
		rotor->T_body_to_hub[0][0] = m_Rotors[ir].Rh(1, 1);
		rotor->T_body_to_hub[1][0] = m_Rotors[ir].Rh(1, 2);
		rotor->T_body_to_hub[2][0] = m_Rotors[ir].Rh(1, 3);
		rotor->T_body_to_hub[0][1] = m_Rotors[ir].Rh(2, 1);
		rotor->T_body_to_hub[1][1] = m_Rotors[ir].Rh(2, 2);
		rotor->T_body_to_hub[2][1] = m_Rotors[ir].Rh(2, 3);
		rotor->T_body_to_hub[0][2] = m_Rotors[ir].Rh(3, 1);
		rotor->T_body_to_hub[1][2] = m_Rotors[ir].Rh(3, 2);
		rotor->T_body_to_hub[2][2] = m_Rotors[ir].Rh(3, 3);
		// orientation of shaft frame in inertial frame
		Mat3x3 Rac(pCraft->GetRCurr()*m_Rh*m_Rotors[ir].Rh);
		// relative rotation between hub and shaft
		Mat3x3 Rhb(Rac.MulTM(m_Rotors[ir].pHub->GetRCurr()));
		Vec3 Psi(RotManip::VecRot(Rhb));
		rotor->azimuthal_offset = 0.;
		rotor->azimuth = Psi(3)*rotor->rotation_dir;
		silent_cout("ModuleCHARM(" << uLabel << "): "
			"rotor " << ir << "/" << num_rotors
			<< " azimuth=" << rotor->azimuth
			<< std::endl);
		doublereal dPsi = 2.*M_PI/rotor->number_blades;
#if 0
		std::cerr << "*** dPsi=" << dPsi << std::endl;
#endif

		for (int ib = 0; ib < rotor->number_blades; ib++) {
			// FIXME: azimuth of blades increases
			m_Rotors[ir].Blades[ib].Rh = m_Rotors[ir].Rh*RotManip::Rot(Vec3(0., 0., rotor->rotation_dir*ib*dPsi))*Rb2h;
#if 0
			std::cerr << "*** m_Rotors[" << ir << "].Blades[" << ib << "].Rh=" << m_Rotors[ir].Blades[ib].Rh << std::endl;
#endif
		}
		rotor->omega100 = m_Rotors[ir].omega100;
#if 0
		rotor->nominal_thrust_coeff = 0.;
#endif
	}

	for (PD::iterator m_data_iter = m_data.begin(); m_data_iter != m_data.end(); m_data_iter++) {
		if (m_data_iter->pRB) {
			int ir = m_data_iter->pRB->iRotor;
			int ib = m_data_iter->pRB->iBlade;
			int ip = m_data_iter->pRB->iFirst + m_data_iter->iOff;

			m_data_iter->pos = &m_wpaircraft.rotors[ir].blades[ib].control_pts[ip][0];
			m_data_iter->vel = &m_wpaircraft.rotors[ir].blades[ib].cp_velocity[ip][0];

			m_data_iter->tangential_velocity_p = &m_wpaircraft.rotors[ir].blades[ib].tangential_velocity[ip];
			m_data_iter->spanwise_lift_p = &m_wpaircraft.rotors[ir].blades[ib].spanwise_lift[ip];

		} else {
			unsigned idx = m_eval_pts.size();
			m_eval_pts.resize(idx + 1);
			
			m_eval_pts[idx].ref_frame = 0;

			m_data_iter->iOff = idx;
			m_data_iter->pos = &m_eval_pts[idx].xyz[0];
			m_data_iter->vel = &m_eval_pts[idx].uvw[0];
		}

#if 1
		// initialize, otherwise valgrind may complain
		m_data_iter->vel[0] = 0.;
		m_data_iter->vel[1] = 0.;
		m_data_iter->vel[2] = 0.;
#endif
	}

	Set_int();

	// TODO: generate "charm.inp"
	{
		struct stat st;
		if (stat("charm.inp", &st) != 0) {
			int save_errno = errno;
			silent_cerr("ModuleCHARM::Init_int: unable to stat file \"charm.inp\" (" << save_errno << ": " << strerror(save_errno) << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	m_chglobal.current_time = pDM->pGetDrvHdl()->dGetTime();
	m_chglobal.time_increment = pDM->pGetDrvHdl()->dGetTimeStep();

	silent_cout("ModuleCHARM::Init_int: current_time=" << m_chglobal.current_time << " time_increment=" << m_chglobal.time_increment << std::endl);

	silent_cout("ModuleCHARM::Init_int: calling initWPModule()" <<std::endl);

	// FIXME: one aircraft only
	initWPModule(1,
		&m_wpaircraft,
		m_eval_pts.size(), m_eval_pts.size(), &m_eval_pts[0],
		1,		// FIXME: areEvalPointsFixed?
		&m_chglobal, &m_wpoptions);

	silent_cout("ModuleCHARM::Init_int: called initWPModule()" <<std::endl);
}

void
ModuleCHARM::Set_int(void)
{
	const Vec3& Xac(pCraft->GetXCurr());
	// p, T not used yet
	doublereal rho, c, p, T;
	GetAirProps(Xac, rho, c, p, T);
	Vec3 Vac(pCraft->GetVCurr());
	// Mat3x3 Rac(pCraft->GetRCurr()*m_Rh); // m_R computed by AssRes()
	Vec3 VBac(m_Rac.MulTV(Vac));
	Vec3 WBac(m_Rac.MulTV(pCraft->GetWCurr()));
#if 0
	// FIXME: Euler angles? 123 or 321?  Rac or Rac^T?
	Vec3 EBac = MatR2EulerAngles123(m_Rac);
#endif

	// update global data
	Vec3 Vinf(0.);
	if (fGetAirVelocity(Vinf, Xac)) {
#if 0
		// Apparently, inertial_winds is incompatible 
		// with INFLOW == 1 (analysis in the aircraft frame)
		m_chglobal.inertial_winds[0] = Vinf(1);
		m_chglobal.inertial_winds[1] = Vinf(2);
		m_chglobal.inertial_winds[2] = Vinf(3);
#endif
		Vac -= Vinf;
	}

	doublereal dTime = pDM->pGetDrvHdl()->dGetTime();
	m_chglobal.current_time = dTime;
	m_chglobal.time_increment = pDM->pGetDrvHdl()->dGetTimeStep();

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
#if 0
	m_wpaircraft.euler_angles[0] = EBac(1);
	m_wpaircraft.euler_angles[1] = EBac(2);
	m_wpaircraft.euler_angles[2] = EBac(3);
#endif
	// T_inertial_to_body = m_Rac^T
	// NOTE: C indexes need to be exchanged
	m_wpaircraft.T_inertial_to_body[0][0] = m_Rac(1, 1);
	m_wpaircraft.T_inertial_to_body[1][0] = m_Rac(1, 2);
	m_wpaircraft.T_inertial_to_body[2][0] = m_Rac(1, 3);
	m_wpaircraft.T_inertial_to_body[0][1] = m_Rac(2, 1);
	m_wpaircraft.T_inertial_to_body[1][1] = m_Rac(2, 2);
	m_wpaircraft.T_inertial_to_body[2][1] = m_Rac(2, 3);
	m_wpaircraft.T_inertial_to_body[0][2] = m_Rac(3, 1);
	m_wpaircraft.T_inertial_to_body[1][2] = m_Rac(3, 2);
	m_wpaircraft.T_inertial_to_body[2][2] = m_Rac(3, 3);

	// TODO: rotor(s)
	for (unsigned ir = 0; ir < m_Rotors.size(); ir++) {
		wpRotorSurface *rotor = &m_wpaircraft.rotors[ir];

		// hub position
		Vec3 Xhb(m_Rac.MulTV(m_Rotors[ir].pHub->GetXCurr() - Xac));
		rotor->hub_position[0] = Xhb(1);
		rotor->hub_position[1] = Xhb(2);
		rotor->hub_position[2] = Xhb(3);

		// azimuth
		// Mat3x3 Rac(pCraft->GetRCurr()*m_Rh*m_Rotors[ir].Rh);	// m_R computed by AssRes()
		Mat3x3 Rac(m_Rac*m_Rotors[ir].Rh);
		Mat3x3 Rhb(m_Rac.MulTM(m_Rotors[ir].pHub->GetRCurr()));
		Vec3 Psi(RotManip::VecRot(Rhb));
		rotor->azimuth = Psi(3);
		while (rotor->azimuth < 0.) {
			rotor->azimuth += 2.*M_PI;
		}
		while (rotor->azimuth > 2.*M_PI) {
			rotor->azimuth -= 2.*M_PI;
		}
		rotor->azimuth = rotor->rotation_dir*Psi(3);

		// rotor velocity
		Vec3 ehb3(m_Rotors[ir].pHub->GetRCurr().GetVec(3));
		rotor->rotor_speed = std::abs(ehb3*(m_Rotors[ir].pHub->GetWCurr() - pCraft->GetWCurr()));

#if 0
		silent_cout("ModuleCHARM(" << GetLabel() << "): "
			"ir=" << ir << " psi=" << rotor->azimuth
			<< " omega=" << rotor->rotor_speed << std::endl);
#endif

#if 0
		for (unsigned ib = 0; ib < m_Rotors[ir].Blades.size(); ib++) {
			wpBlade *blade = &rotor->blades[ib];

			for (unsigned ie = 0; ie < m_Rotors[ir].Blades[ib].Elems.size(); ie++) {
				int iFirst = m_Rotors[ir].Blades[ib].Elems[ie].pRB->iFirst;
				int iCount = m_Rotors[ir].Blades[ib].Elems[ie].pRB->iCount;
				int iIdx = m_Rotors[ir].Blades[ib].Elems[ie].pRB->iIdx;

				ASSERT(iIdx != -1);

				for (int ip = 0; ip < iCount; ip++) {
					blade->tangential_velocity[iFirst + ip] = m_data[iIdx].tangential_velocity;
					blade->spanwise_lift[iFirst + ip] = m_data[iIdx].spanwise_lift;
					blade->control_pts[iFirst + ip][0] = m_data[iIdx].X(1);
					blade->control_pts[iFirst + ip][1] = m_data[iIdx].X(2);
					blade->control_pts[iFirst + ip][2] = m_data[iIdx].X(3);
				}
			}
		}
#endif
	}

#if 0
	std::cerr << "### ModuleCHARM::Update_int: m_data.size=" << m_data.size() << std::endl;
#endif

	// loop on all data points, regardless of their type
	// to set X; blade points also set tangential velocity and lift
	for (PD::iterator i = m_data.begin(); i != m_data.end(); i++) {
		i->X.PutTo(i->pos);
		if (i->pRB) {
			*(i->tangential_velocity_p) = i->tangential_velocity;
			*(i->spanwise_lift_p) = i->spanwise_lift;
		}

#if 0
		if (i->pRB) {
			std::cerr << "Rotors[" << i->pRB->iRotor << "].Blades[" << i->pRB->iBlade << "] X=" << i->X << std::endl;
		} else {
			std::cerr << "X=" << i->X << std::endl;
		}
#endif
	}
}

void
ModuleCHARM::Update_int(void)
{
	Set_int();

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
ModuleCHARM::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	Vec3 V;

#if 0
	std::cerr << "ModuleCHARM(" << GetLabel() << ")::GetInducedVelocity: "
		<< psElemNames[type] << "(" << uLabel << "):" << uPnt << ": "
		"X={" << X << "} iFirstAssembly=" << iFirstAssembly << std::endl;
#endif

	if (iFirstAssembly == 1) {
		bool bGotIt(false);
		for (PD::const_iterator i = m_data.begin(); i != m_data.end(); i++) {
			if (i->type == type && i->label == uLabel && i->counter == uPnt) {

#if 0
				std::cerr << "ModuleCHARM(" << GetLabel() << ")::GetInducedVelocity: "
					<< psElemNames[type] << "(" << uLabel << "):" << uPnt << ": "
					"gotit" << std::endl;
#endif

				bGotIt = true;
				break;
			}
		}

		if (!bGotIt) {
			unsigned idx = m_data.size();
			m_data.resize(idx + 1);
			m_data[idx].type = type;
			m_data[idx].label = uLabel;
			m_data[idx].counter = uPnt;

			// Mat3x3 Rac(pCraft->GetRCurr()*m_Rh);
			m_data[idx].X = m_Rac.MulTV(X - pCraft->GetXCurr());
		}

		V = Zero3;

	} else {
		// NOTE: we only need the final value of the position
		// to update it when the new step of CHARM WP is needed
		// for this purpose it may be preferable to only save X
		// and manipulate it later
		// there is no need to cache the value of the induced
		// velocity, as this would only save a Vec3 constructor
		// from double[3]
		if (m_data_vel_iter->pRB) {
			int ir = m_data_vel_iter->pRB->iRotor;
			int ib = m_data_vel_iter->pRB->iBlade;
			int iRotationDir = m_Rotors[ir].rotation_dir;

			Mat3x3 RTmp(m_Rotors[ir].pHub->GetRCurr()*m_Rotors[ir].Blades[ib].Rh);

#if 0
			std::cerr << "*** ir=" << ir << std::endl;
			std::cerr << "*** ib=" << ib << std::endl;
			std::cerr << "*** m_Rac=" << m_Rac << std::endl;
			std::cerr << "*** m_Rotors[" << ir << "].Blades[" << ib << "].Rh=" << m_Rotors[ir].Blades[ib].Rh << std::endl;
			std::cerr << "*** X={" << RTmp.MulTV(X - m_Rotors[ir].pHub->GetX()) << "}" << std::endl;
#endif

			// blade points (XSIM) are in the blade frame
			// the blade frame originates from root cutout
			// in the disk plane
			// axis 1 points outwards radially
			// axis 2 points from leading to trailing edge
			// axis 3 is 1 cross 2 (down for counter-clockwise)
			m_data_frc_iter->X = RTmp.MulTV(X - m_Rotors[ir].pHub->GetX());
			m_data_frc_iter->X(1) -= m_Rotors[ir].root_cutout;
			m_data_frc_iter->X(2) *= iRotationDir;

			V = RTmp*Vec3(m_data_vel_iter->vel[0], iRotationDir*m_data_vel_iter->vel[1], m_data_vel_iter->vel[2]);

#if 0
			std::cerr << "*** V=" << V << std::endl;
#endif

		} else {
			// set for later...
			m_data_vel_iter->X = m_Rac.MulTV(X - pCraft->GetXCurr());

			V = m_Rac*Vec3(m_data_vel_iter->vel);

#if 0
			std::cerr << "+++ V=" << V << std::endl;
#endif
		}

#if 0
		std::cerr << "ModuleCHARM(" << GetLabel() << ")::GetInducedVelocity: "
			"X={" << X << "} V={" << V << "} Xloc={" << m_data_vel_iter->X <<"}" << std::endl;
#endif

		m_data_vel_iter++;
	}

	return V;
}

void
ModuleCHARM::AddSectionalForce(Elem::Type type,
	unsigned int uLabel, unsigned uPnt,
	const Vec3& F, const Vec3& M, doublereal dW,
	const Vec3& X, const Mat3x3& R, const Vec3& V, const Vec3& W)
{
#if 0
	std::cerr << "ModuleCHARM(" << GetLabel() << ")::AddSectionalForce: "
		<< psElemNames[type] << "(" << uLabel << "):" << uPnt << std::endl;
#endif

	if (iFirstAssembly == 1) {
		unsigned idx = m_data.size();
		ASSERT(idx > 0);
		ASSERT(uLabel != unsigned(-1));
		idx--;

		if (m_data[idx].type != type
			|| m_data[idx].label != uLabel
			|| m_data[idx].counter != uPnt)
		{
			silent_cerr("ModuleCHARM(" << GetLabel() << ")::AddSectionalForce: "
				<< psElemNames[type] << "(" << uLabel << "):" << uPnt << ": "
				"mismatch, expecting "
				<< psElemNames[m_data[idx].type] << "(" << m_data[idx].label << "):" << m_data[idx].counter
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

#if 0
		m_data[idx].label = uLabel;
		if (idx == 0 || (idx > 0 && uLabel != m_data[idx - 1].label)) {
			m_data[idx].counter = 0;
		} else {
			m_data[idx].counter = m_data[idx - 1].counter + 1;
		}
#endif

		m_data[idx].spanwise_lift = 0.;
		m_data[idx].tangential_velocity = 0.;

#if 0
		std::cerr << "ModuleCHARM(" << GetLabel() << ")::AddSectionalForce: "
			<< psElemNames[type] << "(" << uLabel << "):" << uPnt << ": "
			"idx=" << idx << std::endl;
#endif

		// sanity check
		Vec3 Xloc(m_Rac.MulTV(X - pCraft->GetXCurr()));
		if (!m_data[idx].X.IsSame(Xloc, 1e-15)) {
			silent_cerr("ModuleCHARM(" << GetLabel() << ")::AddSectionalForce: "
				<< psElemNames[type] << "(" << uLabel << "):" << uPnt << ": "
				"Xloc = {" << Xloc << "}, GetInducedVelocity({" << m_data[idx].X << "}) "
				"possible mismatch" << std::endl);
		}

		int ir = m_e2b_map[uLabel]->iRotor;
		int ib = m_e2b_map[uLabel]->iBlade;
		int ie = m_e2b_map[uLabel]->iElem;
		if (m_Rotors[ir].Blades[ib].Elems[ie].uLabel != uLabel) {
			silent_cerr("ModuleCHARM(" << GetLabel() << ")::AddSectionalForce: "
				<< psElemNames[type] << "(" << uLabel << "):" << uPnt << ": "
				"mismatch" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (m_e2b_map[uLabel]->iCount == -1) {
			m_e2b_map[uLabel]->iIdx = idx;
			m_e2b_map[uLabel]->iCount = 1;

		} else {
			m_e2b_map[uLabel]->iCount++;
		}

		m_data[idx].pRB = m_e2b_map[uLabel];
		m_data[idx].iOff = m_e2b_map[uLabel]->iCount - 1;

#if 0
		std::cerr << "ModuleCHARM(" << GetLabel() << ")::AddSectionalForce: "
				<< psElemNames[type] << "(" << uLabel << "):" << uPnt << ": "
			<< " ir=" << ir << " ib=" << ib << " ie=" << ie << " num_points=" << m_e2b_map[uLabel]->iCount
			<< std::endl;
#endif

	} else {
		// TODO: maybe this can be computed once for all at each AssRes?
		// Mat3x3 Rac(pCraft->GetRCurr()*m_Rh);

		while (m_data_frc_iter->pRB == 0) {
			if (m_data_frc_iter == m_data.end()) {
				silent_cerr("ModuleCHARM(" << GetLabel() << ")::AddSectionalForce: "
					<< psElemNames[type] << "(" << uLabel << "):" << uPnt << ": "
					"unexpected end of iterator" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			m_data_frc_iter++;
		}

		// resolve force, moment and point in craft's reference frame
		Vec3 Vloc(R.MulTV(V));

#if 0
		std::cerr << "*** Vloc={" << Vloc << "}" << std::endl;
#endif

		Vloc(3) = 0.;

		Vec3 Floc(R.MulTV(F));

#if 0
		std::cerr << "*** Floc={" << Floc << "}" << std::endl;
#endif

		Floc(3) = 0.;
		m_data_frc_iter->tangential_velocity = Vloc.Norm();
		if (m_data_frc_iter->tangential_velocity > std::numeric_limits<doublereal>::epsilon()) {
			Vloc /= m_data_frc_iter->tangential_velocity;
			m_data_frc_iter->spanwise_lift = (Vloc.Cross(Floc))(3);

		} else {
			m_data_frc_iter->spanwise_lift = 0.;
		}

#if 1
		std::cerr << "    dV=" << m_data_frc_iter->tangential_velocity
			<< " dF=" << m_data_frc_iter->spanwise_lift
			<< " X={" << m_data_frc_iter->X << "}" << std::endl;
#endif

		m_data_frc_iter++;
	}
}

void
ModuleCHARM::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	NO_OP;
	// Update_int();
}

void
ModuleCHARM::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP)
{
	// NO_OP;
	if (iFirstAssembly) {
		return;
	}

	Update_int();
}

void
ModuleCHARM::Output(OutputHandler& OH) const
{
#if 0
	std::cerr << "ModuleCHARM(" << GetLabel() << ")::Output: size=" << m_data.size() << std::endl;
#endif

	// should do something useful
	if (fToBeOutput()) {
		std::ostream& out = OH.Loadable();

		for (PD::const_iterator i = m_data.begin(); i != m_data.end(); i++) {
			out << GetLabel() << "#" << i->label << "#" << i->counter
				<< " " << i->spanwise_lift
				<< " " << i->tangential_velocity
				<< " " << i->X
				<< " " << i->vel[0] << " " << i->vel[1] << " " << i->vel[2]
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
#if 0
	std::cerr << "ModuleCHARM(" << GetLabel() << ")::AssJac: iFirstAssembly=" << iFirstAssembly << std::endl;
#endif

	// should do something useful
	WorkMat.SetNullMatrix();

	// re-initialize iterators to loop over point data
	// when AddSectionalForces() and GetInducedVelocity() are called
	m_data_frc_iter = m_data.begin();
	while (m_data_frc_iter != m_data.end() && m_data_frc_iter->pRB == 0) {
		m_data_frc_iter++;
	}

	m_data_vel_iter = m_data.begin();

	return WorkMat;
}

SubVectorHandler& 
ModuleCHARM::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
#if 0
	std::cerr << "ModuleCHARM(" << GetLabel() << ")::AssRes: iFirstAssembly=" << iFirstAssembly << std::endl;
#endif

	m_Rac = pCraft->GetRCurr()*m_Rh;

	if (iFirstAssembly) {
		if (iFirstAssembly == 1) {
			// all topological information should be available;
			// module can be initialized
			Init_int();
		}

		iFirstAssembly--;
	}

	// re-initialize iterators to loop over point data
	// when AddSectionalForces() and GetInducedVelocity() are called
	m_data_frc_iter = m_data.begin();
	while (m_data_frc_iter != m_data.end() && m_data_frc_iter->pRB == 0) {
		m_data_frc_iter++;
	}

	m_data_vel_iter = m_data.begin();

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
	ASSERT(pDM == this->pDM);
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

bool
mbcharm_set(void)
{
	UserDefinedElemRead *rf = new UDERead<ModuleCHARM>;

	bool b = SetUDE("charm", rf);
	if (!b) {
		delete rf;
	}

	return b;
}

// #ifdef STATIC_MODULES, the function is registered by InitUDE()
#ifndef STATIC_MODULES
extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!mbcharm_set()) {
		silent_cerr("module-charm: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}
#endif // ! STATIC_MODULES
