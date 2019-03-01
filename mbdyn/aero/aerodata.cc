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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "mynewmem.h"
#include "aerodata_impl.h"
#include "gauss.h"
#include "submat.h"
#include "aerod2.h"
#include "dataman.h"
#include "drive_.h"

/* AeroMemory - begin */

AeroMemory::AeroMemory(DriveCaller *pt)
: a(0), t(0), iPoints(0), numUpdates(0), pTime(pt)
{
	NO_OP;
}

AeroMemory::~AeroMemory(void)
{
	if (iPoints > 0) {
		if (pTime) {
			SAFEDELETE(pTime);
		}

		if (a) {
			SAFEDELETEARR(a);
		}
	}
}

#define USE_POLCOE
void
AeroMemory::Predict(int i, doublereal alpha, doublereal &alf1, doublereal &alf2)
{
	/* FIXME: this should be s, but I don't want a malloc here */
	doublereal coe[3];
	int s = StorageSize();
#ifdef USE_POLCOE
	/*
	 * FIXME: actually this is not the order of the polynomial,
	 * but the number of coefficients, e.g. a second-order polynomial
	 * has three coefficients :)
	 */
	integer		order = s;
#endif /* USE_POLCOE */

	ASSERT(s > 0);
	ASSERT(i < iPoints);

	doublereal *aa = a + s*i;
	doublereal *tt = t + s*i;

	aa[s-1] = alpha;
	tt[s-1] = pTime->dGet();

	if (numUpdates >= s) {
#ifdef USE_POLCOE
		__FC_DECL__(polcoe)(tt, aa, &order, coe);
#else /* !USE_POLCOE */
		coe[2] = ((aa[2]-aa[0])/(tt[2]-tt[0])
				- (aa[1]-aa[0])/(tt[1]-tt[0]))/(tt[2]-tt[1]);
		coe[1] = (aa[1]-aa[0])/(tt[1]-tt[0]) - (tt[1]+tt[0])*coe[2];

#if 0	/* not used ... */
		coe[1] = aa[0] - tt[0]*coe[1] - tt[0]*tt[0]*coe[2];
#endif
#endif /* !USE_POLCOE */

#if 0
			std::cerr << "aa[0:2]= " << aa[0] << "," << aa[1] << "," << aa[2] << std::endl
				<< "tt[0:2]= " << tt[0] << "," << tt[1] << "," << tt[2] << std::endl
				<< "coe[0:2]=" << coe[0] << "," << coe[1] << "," << coe[2] << std::endl;
#endif /* 0 */

		alf1 = coe[1]+2.*coe[2]*tt[2];
		alf2 = 2.*coe[2];
	} else {
		alf1 = 0.;
		alf2 = 0.;
	}
}

void
AeroMemory::Update(int i)
{
	int s = StorageSize();
	if (s > 0) {
		/*
		 * shift back angle of attack and time
		 * for future interpolation
		 */
		doublereal *aa = a + s*i;
		doublereal *tt = t + s*i;

		ASSERT(i < iPoints);

		for (int j = 1; j < s; j++) {
			aa[j-1] = aa[j];
			tt[j-1] = tt[j];
		}

		if (i == 0) {
			numUpdates++;
		}
	}
}

void
AeroMemory::SetNumPoints(int i)
{
	int s = StorageSize();

	iPoints = i;

	if (s > 0) {
		SAFENEWARR(a, doublereal, 2*s*iPoints);
		t = a + s*iPoints;

#ifdef HAVE_MEMSET
		memset(a, 0, 2*s*iPoints*sizeof(doublereal));
#else /* !HAVE_MEMSET */
		for (int i = 0; i < 2*s*iPoints; i++) {
			a[i] = 0.;
		}
#endif /* !HAVE_MEMSET */
	}
}

int
AeroMemory::GetNumPoints(void) const
{
	return iPoints;
}

/* AeroMemory - end */

/* C81Data - begin */

C81Data::C81Data(unsigned int uLabel)
: WithLabel(uLabel)
{
	NO_OP;
}

/* C81Data - end */

/* AeroData - begin */

AeroData::AeroData(int i_p, int i_dim,
	AeroData::UnsteadyModel u, DriveCaller *ptime)
: AeroMemory(ptime), unsteadyflag(u), Omega(0.)
{
	// silence static analyzers
	VAM.density = -1.;
	VAM.sound_celerity = -1.;
	VAM.chord = -1.;
	VAM.force_position = 0.;
	VAM.bc_position = 0.;
	VAM.twist = 0.;

	if (u != AeroData::STEADY) {
		if (ptime == 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	SetNumPoints(i_p*i_dim);
}

AeroData::~AeroData(void)
{
	NO_OP;
}

AeroData::UnsteadyModel
AeroData::Unsteady(void) const
{
	return unsteadyflag;
}

void
AeroData::SetAirData(const doublereal& rho, const doublereal& c)
{
	VAM.density = rho;
	VAM.sound_celerity = c;
}

int
AeroData::StorageSize(void) const
{
	switch (unsteadyflag) {
	case AeroData::HARRIS:
	case AeroData::BIELAWA:
		return 3;

	case AeroData::STEADY:
		return 0;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
AeroData::SetSectionData(const doublereal& abscissa,
	const doublereal& chord,
	const doublereal& forcepoint,
	const doublereal& velocitypoint,
	const doublereal& twist,
	const doublereal& omega)
{
	VAM.chord = chord;
	VAM.force_position = forcepoint;
	VAM.bc_position = velocitypoint;
	VAM.twist = twist;
	Omega = omega;
}

std::ostream&
AeroData::RestartUnsteady(std::ostream& out) const
{
	switch (unsteadyflag) {
	case AeroData::HARRIS:
		out << ", unsteady, harris";
		break;

	case AeroData::BIELAWA:
		out << ", unsteady, bielawa";
		break;

	default:
		break;
	}

	return out;
}

int
AeroData::GetForcesJacForwardDiff_int(int i, const doublereal* W, doublereal* TNG0, Mat6x6& J, outa_t& OUTA)
{
	const doublereal epsilon = 1.e-3;
	const doublereal nu = 1.e-9;

	Vec3 V(&W[0]);
	Vec3 Omega(&W[3]);

	doublereal dv = V.Norm();
	doublereal dw = Omega.Norm();

	doublereal TNG[6];

	GetForces(i, W, TNG0, OUTA);

	doublereal *WW = const_cast<doublereal *>(W);
	for (unsigned int iColm1 = 0; iColm1 < 6; iColm1++)	{
		doublereal dorig;
		doublereal delta;

		dorig = WW[iColm1];
		if (iColm1 < 3) {
			delta = dv*epsilon + nu;

		} else {
			delta = dw*epsilon + nu;
		}
		WW[iColm1] = dorig + delta;

		GetForces(i, WW, TNG, OUTA);

		for (unsigned int iRowm1 = 0; iRowm1 < 6; iRowm1++) {
			J.Put(iRowm1 + 1, iColm1 + 1, (TNG[iRowm1] - TNG0[iRowm1])/delta);
		}

		WW[iColm1] = dorig;
	}

	return 0;
}

int
AeroData::GetForcesJacCenteredDiff_int(int i, const doublereal* W, doublereal* TNG0, Mat6x6& J, outa_t& OUTA)
{
	const doublereal epsilon = 1.e-3;
	const doublereal nu = 1.e-9;

	Vec3 V(&W[0]);
	Vec3 Omega(&W[3]);

	doublereal dv = V.Norm();
	doublereal dw = Omega.Norm();

	doublereal TNGp[6];
	doublereal TNGm[6];

	GetForces(i, W, TNG0, OUTA);

	doublereal *WW = const_cast<doublereal *>(W);
	for (unsigned int iColm1 = 0; iColm1 < 6; iColm1++)	{
		doublereal dorig;
		doublereal delta;

		dorig = WW[iColm1];
		if (iColm1 < 3) {
			delta = dv*epsilon + nu;

		} else {
			delta = dw*epsilon + nu;
		}

		WW[iColm1] = dorig + delta;

		GetForces(i, WW, TNGp, OUTA);

		WW[iColm1] = dorig - delta;

		GetForces(i, WW, TNGm, OUTA);

		doublereal delta2 = 2.*delta;
		for (unsigned int iRowm1 = 0; iRowm1 < 6; iRowm1++) {
			J.Put(iRowm1 + 1, iColm1 + 1, (TNGp[iRowm1] - TNGm[iRowm1])/delta2);
		}

		WW[iColm1] = dorig;
	}

	return 0;
}

unsigned int
AeroData::iGetNumDof(void) const
{
	return 0;
}

DofOrder::Order
AeroData::GetDofType(unsigned int) const
{
	silent_cerr("AeroData: "
		"GetDofType() is undefined because aerodynamic data "
		"has no degrees of freedom" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

int
AeroData::GetForces(int i, const doublereal* W, doublereal* TNG, outa_t& OUTA)
{
	silent_cerr("AeroData: GetForces() is undefined" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

int
AeroData::GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
	silent_cerr("AeroData: GetForcesJac() is undefined" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
AeroData::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr,
	integer iFirstIndex, integer iFirstSubIndex,
	int i, const doublereal* W, doublereal* TNG, outa_t& OUTA)
{
	silent_cerr("AeroData: AssRes() is undefined" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
AeroData::AssJac(FullSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr,
	integer iFirstIndex, integer iFirstSubIndex,
	const Mat3xN& vx, const Mat3xN& wx, Mat3xN& fq, Mat3xN& cq,
	int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
	silent_cerr("AeroData: AssJac() is undefined" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

void
AeroData::AfterConvergence(
	int i,
	const VectorHandler& X,
	const VectorHandler& XP)
{
	NO_OP;
}

/* AeroData - end */

static AeroData::UnsteadyModel
ReadUnsteadyFlag(MBDynParser& HP)
{
	if (HP.IsArg()) {
		AeroData::UnsteadyModel eInst = AeroData::STEADY;
		if (HP.IsKeyWord("unsteady")) {
			/*
			 * swallow "unsteady" keyword
			 */
			if (HP.IsKeyWord("steady")) {
				eInst = AeroData::STEADY;
			} else if (HP.IsKeyWord("harris")) {
				eInst = AeroData::HARRIS;
			} else if (HP.IsKeyWord("bielawa")) {
				eInst = AeroData::BIELAWA;
			} else {
				/* demote to pedantic, because the integer
				 * form allows to change unsteady model
				 * parametrically (while waiting for string
				 * vars) */
				pedantic_cerr("deprecated unsteady model "
					"given by integer number;"
					" use \"steady\", \"Harris\" or \"Bielawa\" "
					"instead, at line " << HP.GetLineData()
					<< std::endl);

				int i = HP.GetInt();
				if (i < AeroData::STEADY || i >= AeroData::LAST) {
					silent_cerr("illegal unsteady flag "
							"numeric value " << i
							<< " at line "
							<< HP.GetLineData()
							<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				eInst = AeroData::UnsteadyModel(i);
			}
		}

		switch (eInst) {
		case AeroData::STEADY:
		case AeroData::BIELAWA:
			break;

		case AeroData::HARRIS:
			silent_cerr("\"Harris\" unsteady aerodynamics "
					"are not available at line "
					<< HP.GetLineData()
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		default:
	 		silent_cerr("illegal unsteady flag at line "
					<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/*
		 * unsteady flag
		 */
		return eInst;
	}

	/*
	 * default: no unsteady ...
	 */
	return AeroData::STEADY;
}

static void
ReadC81MultipleAeroData(DataManager* pDM, MBDynParser& HP, AeroData** aerodata,
	int iGP, int iDim, bool bInterp = false)
{
	doublereal dcltol = 1.e-6;
	if (bInterp && HP.IsKeyWord("tolerance")) {
		dcltol = HP.GetReal();
		if (dcltol <= 0.) {
			silent_cerr("ReadC81MultipleAeroData: "
				"invalid c81 data tolerance at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	integer nProfiles = HP.GetInt();
	if (nProfiles <= 0) {
		silent_cerr("Need at least one airfoil at line "
				<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	std::vector<unsigned> profiles(nProfiles);
	std::vector<doublereal> upper_bounds(nProfiles);
	std::vector<const c81_data *> data(nProfiles);

	for (int i = 0; i < nProfiles; i++) {
		profiles[i] = (unsigned)HP.GetInt();
		upper_bounds[i] = HP.GetReal();
		if (upper_bounds[i] <= -1.) {
			if (!bInterp || upper_bounds[i] < -1.) {
				silent_cerr("upper bound "
					<< i + 1 << " = " << upper_bounds[i]
					<< " too small at line "
					<< HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

		} else if (upper_bounds[i] > 1.) {
			silent_cerr("upper bound "
				<< i + 1 << " = "
				<< upper_bounds[i]
				<< " too large at line "
				<< HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (i > 0 && upper_bounds[i] <= upper_bounds[i-1]) {
			silent_cerr("upper bound "
				<< i+1 << " = " << upper_bounds[i]
				<< " not in increasing order "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		data[i] = HP.GetC81Data(profiles[i]);
		if (data[i] == 0) {
			silent_cerr("Unable to find airfoil "
				<< profiles[i] << " at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
 		DEBUGLCOUT(MYDEBUG_INPUT, "airfoil data " << i+1
			<< " is from file c81 " << profiles[i]
			<< std::endl);
	}

	if (upper_bounds[nProfiles - 1] != 1.) {
		silent_cerr("warning: the last upper bound should be 1.0 "
			"at line " << HP.GetLineData() << std::endl);
	}

	AeroData::UnsteadyModel eInst = ReadUnsteadyFlag(HP);
	DriveCaller *ptime = 0;
	if (eInst != AeroData::STEADY) {
		SAFENEWWITHCONSTRUCTOR(ptime,
			TimeDriveCaller,
			TimeDriveCaller(pDM->pGetDrvHdl()));
	}

	if (bInterp) {
		SAFENEWWITHCONSTRUCTOR(*aerodata,
			C81InterpolatedAeroData,
			C81InterpolatedAeroData(iGP, iDim,
				eInst, profiles, upper_bounds, data,
				dcltol, ptime));

	} else {
		SAFENEWWITHCONSTRUCTOR(*aerodata,
			C81MultipleAeroData,
			C81MultipleAeroData(iGP, iDim,
				eInst, profiles, upper_bounds, data,
				ptime));
	}
}

void
ReadAeroData(DataManager* pDM, MBDynParser& HP, int iDim,
	Shape** ppChord, Shape** ppForce,
	Shape** ppVelocity, Shape** ppTwist, Shape** ppTipLoss,
	integer* piNumber, DriveCaller** ppDC,
	AeroData** aerodata)
{
	DEBUGCOUTFNAME("ReadAeroData");

	/* Keywords */
	const char* sKeyWords[] = {
		"naca0012",
		"rae9671",
		"c81",

		"theodorsen",

		0
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,
		NACA0012 = 0,
		RAE9671,
		C81,

		THEODORSEN,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	*ppChord = ReadShape(HP);
	*ppForce = ReadShape(HP);
	*ppVelocity = ReadShape(HP);
	*ppTwist = ReadShape(HP);

	*ppTipLoss = 0;
	if (HP.IsKeyWord("tip" "loss")) {
		*ppTipLoss = ReadShape(HP);

	} else {
		SAFENEWWITHCONSTRUCTOR(*ppTipLoss, ConstShape1D,
			ConstShape1D(1.));
	}

	*piNumber = HP.GetInt();
	if (*piNumber <= 0) {
		silent_cerr("need at least 1 Gauss integration point at line "
				<< HP.GetLineData()  << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	DEBUGLCOUT(MYDEBUG_INPUT, "Gauss points number: "
			<< *piNumber << std::endl);

	if (HP.IsKeyWord("control")) {
		/* Driver di un'eventuale controllo */
		*ppDC = HP.GetDriveCaller();

	} else {
		SAFENEW(*ppDC, NullDriveCaller);
	}

	if (HP.IsArg()) {
		switch (HP.IsKeyWord()) {
 		default:
	  		silent_cerr("unknown airfoil type \"" << HP.GetString()
					<< "\" at line " << HP.GetLineData()
					<< "; using default (NACA0012)"
					<< std::endl);

 		case NACA0012: {
	  		DEBUGLCOUT(MYDEBUG_INPUT,
				   "airfoil is NACA0012" << std::endl);

			AeroData::UnsteadyModel eInst = ReadUnsteadyFlag(HP);
			DriveCaller *ptime = 0;
			if (eInst != AeroData::STEADY) {
				SAFENEWWITHCONSTRUCTOR(ptime, TimeDriveCaller,
						TimeDriveCaller(pDM->pGetDrvHdl()));
			}
	  		SAFENEWWITHCONSTRUCTOR(*aerodata,
				STAHRAeroData,
				STAHRAeroData(*piNumber, iDim, eInst, 1, ptime));
	  		break;
 		}

 		case RAE9671: {
	  		DEBUGLCOUT(MYDEBUG_INPUT,
				"airfoil is RAE9671" << std::endl);

			AeroData::UnsteadyModel eInst = ReadUnsteadyFlag(HP);
			DriveCaller *ptime = 0;
			if (eInst != AeroData::STEADY) {
				SAFENEWWITHCONSTRUCTOR(ptime, TimeDriveCaller,
						TimeDriveCaller(pDM->pGetDrvHdl()));
			}
	  		SAFENEWWITHCONSTRUCTOR(*aerodata,
				STAHRAeroData,
				STAHRAeroData(*piNumber, iDim, eInst, 2, ptime));
	  		break;
 		}

	 	/*
		 * uso tabelle standard, che vengono cercate in base all'indice
		 * (sistemare)
		 */
 		case C81:
			if (HP.IsKeyWord("multiple")) {
				ReadC81MultipleAeroData(pDM, HP, aerodata, *piNumber, iDim);

			} else if (HP.IsKeyWord("interpolated")) {
				ReadC81MultipleAeroData(pDM, HP, aerodata, *piNumber, iDim, true);

			} else {
	  			unsigned iProfile = (unsigned)HP.GetInt();
		  		const c81_data* data = HP.GetC81Data(iProfile);

		  		DEBUGLCOUT(MYDEBUG_INPUT,
					"airfoil data is from file c81 "
					<< iProfile << std::endl);
				AeroData::UnsteadyModel
					eInst = ReadUnsteadyFlag(HP);
				DriveCaller *ptime = 0;
				if (eInst != AeroData::STEADY) {
					SAFENEWWITHCONSTRUCTOR(ptime,
							TimeDriveCaller,
							TimeDriveCaller(pDM->pGetDrvHdl()));
				}
	  			SAFENEWWITHCONSTRUCTOR(*aerodata,
					C81AeroData,
					C81AeroData(*piNumber, iDim,
						eInst, iProfile,
						data, ptime));
			}
			break;

		case THEODORSEN: {
			if (!HP.IsKeyWord("c81")) {
				silent_cerr("Theodorsen aerodata: warning, assuming \"c81\" at line "
				<< HP.GetLineData() << std::endl);
			}

			AeroData *pa = 0;
			if (HP.IsKeyWord("multiple")) {
				ReadC81MultipleAeroData(pDM, HP, &pa, *piNumber, iDim);

			} else if (HP.IsKeyWord("interpolated")) {
				ReadC81MultipleAeroData(pDM, HP, &pa, *piNumber, iDim, true);

			} else {
	  			unsigned iProfile = (unsigned)HP.GetInt();
		  		const c81_data* data = HP.GetC81Data(iProfile);

		  		DEBUGLCOUT(MYDEBUG_INPUT,
					"airfoil data is from file c81 "
					<< iProfile << std::endl);
				AeroData::UnsteadyModel
					eInst = ReadUnsteadyFlag(HP);
				DriveCaller *ptime = 0;
				if (eInst != AeroData::STEADY) {
					SAFENEWWITHCONSTRUCTOR(ptime,
							TimeDriveCaller,
							TimeDriveCaller(pDM->pGetDrvHdl()));
				}
	  			SAFENEWWITHCONSTRUCTOR(pa,
					C81AeroData,
					C81AeroData(*piNumber, iDim,
						eInst, iProfile,
						data, ptime));
			}

			DriveCaller *ptime = 0;
			SAFENEWWITHCONSTRUCTOR(ptime,
					TimeDriveCaller,
					TimeDriveCaller(pDM->pGetDrvHdl()));

			SAFENEWWITHCONSTRUCTOR(*aerodata,
				TheodorsenAeroData,
				TheodorsenAeroData(*piNumber, iDim,
					pa, ptime));
			} break;
		}

	} else {
		/* FIXME: better abort! */
		silent_cerr("missing airfoil type at line "
				<< HP.GetLineData()
				<< "; using default (NACA0012)"
				<< std::endl);

		AeroData::UnsteadyModel eInst = ReadUnsteadyFlag(HP);
		SAFENEWWITHCONSTRUCTOR(*aerodata,
			STAHRAeroData, STAHRAeroData(*piNumber, iDim, eInst, 1));
	}
}

