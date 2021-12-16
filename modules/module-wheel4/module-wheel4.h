/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2010
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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
 * Author:	Louis Gagnon <louis.gagnon.10@ulaval.ca>
 *		Departement de genie mecanique
 *		Universite Laval
 *		http://www.gmc.ulaval.ca
 */

#ifndef MODULE_WHEEL4_H
#define MODULE_WHEEL4_H

#include "dataman.h"
#include "userelem.h"
#include "drive_.h" // to get current time (may not be needed in the end)

//#include "joint.h" // to get information on viscoelastic elements

class Wheel4
: virtual public Elem, public UserDefinedElem
{
private:

	// patch degrees of freedom,
	doublereal dXx, dXy, dVx, dVy, dVPx, dVPy, dXPx, dXPy;
	// ,patch degrees of freedom


	// wheel node
	StructNode *pWheel;
	Elem *pWheelB; // body
	Elem *pWheelE; // user loadable elem

	// ring node
	StructNode *pRing;
	Elem *pRingB; //body

	// ring node
	StructNode *pPatch;


	// wheel axle direction (wrt/ wheel node)
	Vec3 WheelAxle;

	// (flat) ground node
	StructNode *pGround;



	// friction data
	bool bSlip;
	bool bLoadedRadius;
	DriveCaller *pMuX0;
	DriveCaller *pMuY0;
	DriveCaller *pTr;
	DriveCaller *pMzr;
//	DriveCaller *pMuY1;

	// centripetal force calculation
	doublereal rRatio; // ratio btwn portion of ring in contact with patch and the total ring
	doublereal Krz; // vertical wheel-ring stiffness

	//time variant coefficients
	DriveCaller *pKpa;
	doublereal dKpa;
	Vec3 Kpatv;
	DriveCaller *pCpa;
	doublereal dCpa;
	Vec3 Cpatv;

	//road elevation
	DriveCaller *pRoad;
	doublereal dRoad;
	doublereal dRoadAhead; // for two-point follower
	doublereal dRoadBehind; // for two-point follower
	doublereal dRoadInitial; // road height for simulation onset
	doublereal dRoadPrev;
	doublereal dRoadPrevPrev;

	// thresholds
	doublereal TRH;		// prevents division by zero at null x-velocity at the price of losing validity for velocities above -TRH*1.1 and below TRH*1.1
	doublereal TRHA;	// buffer used to prevent division by zero
	doublereal TRHT;	// prevents division by zero for computation of the angle of the car (and wheels)
	doublereal TRHTA;	// buffer used on angle zero division prevention
	doublereal TRHC;	// cap on kappa
	doublereal TdLs;	// minimum value (cap) for dLs (half-length of two-point follower)
	doublereal TdReDiv;	// minimum value (cap) for divider used to find R_e (effective rolling radius)
	doublereal TdR_e;	// minimum value (cap) for R_e (effective rolling radius)
	doublereal RDA;		// x-pos prior to which road profile is null and after which is interpolated between 0 and value at RDB
	doublereal RDB;		// x-pos where road starts to be taken as fed by driver but using x=0 @ RDB
	doublereal RDL;		// after x reaches RDB+RDL, the road profile will loop over RDL
	int RDLC;			// number of times to deduct x pos from current pos for looping

	// tire data,
	Vec3 Xpa; // elongation of the contact patch spring
	Vec3 Xpar;  //relative to point on ring
	Vec3 distM; // the real physical distance between center of ring and patch after rotation of the patch to respect the slope
	doublereal ddistM; // the magnitude of the home-made loaded radius distM
	doublereal dEffRad; // wheel effective radius
	Vec3 EffRad; // physical (artificial rotation done) vector distance between wheel center and patch
	doublereal R_e; // wheel effective radius computed only for output
	Vec3 XparPrev;  //relative to point on ring, previous timestep
	Vec3 XpaPrev; // elongation of the contact patch spring previous timestep
	Vec3 XparPrevPrev;  //relative to point on ring, previous previous timestep
	Vec3 XpaPrevPrev; // elongation of the contact patch spring, previous previous timestep
	Vec3 Kpa; // stiffness
	Vec3 Cpa; // damping
	Vec3 XpaBC; // elongation (POSITION?) of contact patch before convergence is confirmed
	Vec3 VpaBC; // elongation velocity of contact patch before convergence is confirmed
	Vec3 RpaBC; // angle of contact patch before convergence is confirmed
	Vec3 WpaBC; // angular of contact patch before convergence is confirmed
	Vec3 Vpa;
	Vec3 Vpar; // relative to point on ring
	Vec3 VparWheel; // relative vel betwn patch and wheel
	Vec3 VpaPrev;
	Vec3 Fint; // viscoelastic element force between ring and patch
	Vec3 Fint_old; // to keep value of Fint prior to rotation back into abs. ref. frame
	Vec3 Fint_ring; // viscoelastic element force between ring and patch used for artificial patch rotation
	Vec3 Fpatch; // force on patch
	Vec3 Mint; // viscoelastic element moment between ring and patch
	doublereal Mpa;
	doublereal tr; // pneumatic trail
	doublereal S_ht; // horizontal shift of pneumatic trail
	doublereal S_hf; // horizontal shift of residual torque (equal horiz shift of lateral force)
	doublereal M_zr; // residual torque
	doublereal dt; //timestep
	bool bSwift;
	doublereal curTime;	//current time
	doublereal oldTime; //time at previous timestep
	bool bFirstAC;	// first timestep after convergence (will not reset if timestep changes)
	bool bFirstAP;	// first iter after prediction (resets if timestep changed)
	DriveOwner	tdc;		// time drive
	Vec3 zZero; // to remove a z component of a vector
	Vec3 pcRing; // contact point of ring
	Vec3 pcRingPrev; // contact point of ring previous timestep
	Vec3 RingRad; // vector radius of ring
	Vec3 RingRadPrev; // previous vector radius of ring
	Vec3 VpcRingPrev;
	Vec3 VpcRing;
	Vec3 VpcRingPrevPrev;
	Vec3 Xring; // abs pos of ring axle
	Vec3 Vwheel; // current velocity of wheel
	Vec3 fwd; // forward direction of wheel in absol. ref frame
	Vec3 fwdRing; // forward direction of ring in absol. ref frame
	Vec3 fwdRingFlat; // forward direction of ring in absol. ref frame disregarding inclination of the slope
	Vec3 lat; // lateral direction of wheel in absol. ref frame
	Vec3 latRing; // lateral direction of ring in absol. ref frame
	Vec3 latRingFlat; // lateral direction of ring in absol. ref frame disregarding inclination of the slope
	Vec3 n; // ground orientation in the absolute frame
	Vec3 nPrev; // ground orientation in the absolute frame previous timestep
	Vec3 fwdRingPrev; // fwd unit vector of the ring node at previous timestep
	doublereal dRoadVel; // road velocity in the z-dir
	Vec3 Xparp; // rotated relative displacement vector of the patch prior to multiplication by the stiffness
	Vec3 Vparp; // rotated relative velocity vector of the patch prior to multiplication by the viscosity
	doublereal Fn; // force on patch in z-direction
	doublereal Fcent; // centrifugal force calculated from tire properties and velocity and "applied to patch"
	doublereal dCt; // calculated displacement induced by centrifugal force (calculated before Fcent)
	Vec3 Fr; // rolling resistance force vector that points in the direction of the ring
	bool boolFn; // null if Fn < 0 to help get a proper jacobian
	Vec3 i; // unit vector in x-dir
	Vec3 j; // unit vector in y-dir
	Vec3 k; // unit vector in z-dir
	doublereal dLs; // half-width of the tandem elliptical cam follower
	doublereal dPls; // ratio of lengths between two point follower and contact patch length
	doublereal dR_a1; // r_a1 contact length parameter from Besselink eq. 4.85
	doublereal dR_a2; // r_a2 contact length parameter from Besselink eq. 4.85
	doublereal dLsProj; // half-legnth in abs x-dir of the two point follower projected on the abs x dir
	doublereal dXxProj; // projected x-pos of patch according to road inclination
	doublereal dXxProjPrev; // previous projected x-pos of patch, used to predict the next position
	doublereal dLsProjPrev; // previous projected half length of contact patch, used to predict the next position
	doublereal dt_maxF; // maximum divison of dt per timestep
	doublereal dt_minF; // minimum division of dt per timestep
	doublereal dLsProjRatio; // ratio of length of dLs when projected, used for timestep size calculation
	doublereal dtPrev; // previous timestep, used to predict the next position
	doublereal dt_maxH; // max height of step in vertical direction wanted, drive the timestep control
	doublereal dt_adjFactor; // // factor by which the current time step is too big for the bump to come in dt_numAhead times the previous step distance
	doublereal dt_fNow; // factor to apply now to current timestep (ie: divide current timestep by this for the next one)
	doublereal dt_maxstep; // maximum timestep given to input file, needed by timestep driver
	doublereal dt_minstep; // minimum timestep given to input file, needed by timestep driver
	bool dt_On; // bool to enable or disable adjustable timestep calculation
	doublereal dt_Res; // resolution of bump search, keep sufficiently smaller than ring radius (ie: look at every dt_Res meters for a bump in front of the tire) NOTE: will not look behind because it is assumed that the model is not made to move rearwards and that a small move rearward would not influence the timestep because it should already have been adjusted for the bump when driving towards it)
	doublereal dtMax; // maximum recommended timestep
	int dt_numAhead; // number of steps by which to look ahead for a bump in the road profile (timestep controller, not defined by user)
	int dt_minStepsCycle; // minimum number of steps wanted in a force cycle (approx., influences dt)
	doublereal dt_divF; // factor by which to divide the timestep if the force oscillates more than wanted (as determined by TminS)
	doublereal dt_divF3; // factor if the force changes 3 times of sign in the last dt_minStepsCycle steps
	doublereal dt_divF4; // factor if the force changes 4 or more times of sign in the last dt_minStepsCycle steps
	std::vector< std::vector<int> > FintSignCk;
	Vec3 FintPrev; // previous Fint for timestep control
	Vec3 FintPrevPrev; // previous-previous Fint for timestep control
	Vec3 XparpPrev; // previous Xpar for timestep control
	Vec3 XparpPrevPrev; // previous-previous Xparp for timestep control
	doublereal dn; // magnitude square of normal vector to road (n)
	doublereal dvx; // relative speed between center of wheel and contact point on tire in the forward direction
	doublereal dvax; // speed of axle in the forward direction
	doublereal dvay; // speed of axle in the lateral direction
	Vec3 va; 	// relative speed between wheel axle and ground
	doublereal dXxPrev; // x-pos of patch at previous timestep
	doublereal E;   // total system kin + pot energy
	doublereal KE;   // total system kin energy
	doublereal PE;   // total system pot energy
	// ,tire data



	// variables used by the Jacobian,

	doublereal derivSign; // ensures the proper derivative sign by respecting the abs(dvax) present function.
	bool latBool; // determines whether lateral forces will contribute to the Jacobian
	bool fwdBool; // determines whether longitudinal forces will contribute to the Jacobian

	// ,variables used by the Jacobian

	// output data
	Vec3 F;
	Vec3 M;
	Vec3 Mz;
	doublereal dR_0; // Rigid ring radius
	doublereal deltaPrev; // tire deflection, from previous timestep
	doublereal dSr; // long. slip ratio
	doublereal dSa; // lateral slip angle
	doublereal dAlpha;
	doublereal dAlpha_t;
	doublereal dAlpha_r;
	doublereal dMuX;
	doublereal dMuY;
	doublereal q_sy1; // should be between 0.01 and 0.02 usually...
	doublereal q_sy3; //
	doublereal dvao; // reference velocity for rolling resistance velocity influence factor

        bool firstRes; // this_is_the_first_residual 

    doublereal dDebug; // onloy used for debugging output


//	secant of a scalar
        doublereal sec(doublereal x) const {
        		return 1.0 / cos(x);
        };

        int sign(const doublereal x) {
        	if (x >= 0.) {
        		return 1;
        	} else if (x < 0.) {
        		return -1;
        	}
        	return 0;
        };

#ifdef USE_NETCDF
	MBDynNcVar Var_Fint;
	MBDynNcVar Var_Xpar;
	MBDynNcVar Var_Xparp;
	MBDynNcVar Var_dXxProj;
	MBDynNcVar Var_dRoad;
	MBDynNcVar Var_F;
	MBDynNcVar ar_Fn;
	MBDynNcVar Var_debug;
	MBDynNcVar Var_dSr;
	MBDynNcVar Var_ddistM;
	MBDynNcVar Var_Fcent;
	MBDynNcVar Var_dLs;
	MBDynNcVar Var_R_e;
	MBDynNcVar Var_dSa;
	MBDynNcVar Var_dvax;
	MBDynNcVar Var_dvx;
	MBDynNcVar Var_dvay;
	MBDynNcVar Var_dMuY;
	MBDynNcVar Var_dMuX;
	MBDynNcVar Var_KE;
	MBDynNcVar Var_PE;
	MBDynNcVar Var_E;
	MBDynNcVar Var_dRoadAhead;
	MBDynNcVar Var_dRoadBehind;
	MBDynNcVar Var_dCt;
	MBDynNcVar Var_M;
	MBDynNcVar Var_distM;
	MBDynNcVar Var_n;
	MBDynNcVar Var_Xpa;
	MBDynNcVar Var_Vpa;
	MBDynNcVar Var_Vpar;
	MBDynNcVar Var_fwd;
	MBDynNcVar Var_fwdRing;
	MBDynNcVar Var_fwdRingFlat;
	MBDynNcVar Var_pcRing;
	MBDynNcVar Var_VparWheel;
	MBDynNcVar Var_Fr;
	MBDynNcVar Var_Mz;

#endif /* USE_NETCDF */

public:
	Wheel4(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~Wheel4(void);

	virtual void OutputPrepare(OutputHandler &OH);
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

//    /* enables access to private data */
	unsigned int iGetNumPrivData(void) const;
	unsigned int iGetPrivDataIdx(const char *s) const;
	doublereal dGetPrivData(unsigned int i) const;

	int iGetNumConnectedNodes(void) const;
	virtual void SetInitialValue(VectorHandler& XCurr);
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);


	virtual unsigned int iGetNumDof(void) const;
	virtual DofOrder::Order GetDofType(unsigned int i) const;
	virtual DofOrder::Order GetEqType(unsigned int i) const;


	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		      const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

//    /* enables access to private data */
//    virtual unsigned int iGetNumPrivData(void) const;
//    virtual unsigned int iGetPrivDataIdx(const char *s) const;
//    virtual doublereal dGetPrivData(unsigned int i) const;

	void AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP); // to call function after convergence of current step
	void AfterPredict(VectorHandler& X, VectorHandler& XP); // at beginning of iteration


    void CalculateR_e();
    doublereal CapLoop(doublereal Xuncapped) const;
//    void NetCDFPrepare(OutputHandler &OH, char &buf) const;


};

#endif // ! MODULE_WHEEL4_H


/* Wheel4 - end */

/* TimeStep - begin */
#ifndef MODULE_TIMESTEP_H
#define MODULE_TIMESTEP_H

class TimeStep
: virtual public Elem, public UserDefinedElem {
private:

	Elem *pWheelE; // user loadable elem
	std::vector<Elem *> pWheelsE;
//	DriveOwner	tdc;		// time drive


public:
	TimeStep(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~TimeStep(void);

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

	//    /* enables access to private data */
		unsigned int iGetNumPrivData(void) const;
		unsigned int iGetPrivDataIdx(const char *s) const;
		doublereal dGetPrivData(unsigned int i) const;

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

#endif // ! MODULE_TIMESTEP_H
