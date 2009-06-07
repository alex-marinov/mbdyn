/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "mynewmem.h"
// #define USE_C81INTERPOLATEDAERODATA
#include "aerodata_impl.h"
#include "gauss.h"
#include "submat.h"
#include "aerod2.h"

/* STAHRAeroData - begin */

STAHRAeroData::STAHRAeroData(AeroData::UnsteadyModel u, integer p,
		DriveCaller *ptime)
: AeroData(u, ptime), profile(p)
{
	ASSERT(u != AeroData::STEADY ? (ptime != 0) : 1);
}

STAHRAeroData::~STAHRAeroData(void)
{
	NO_OP;
}

std::ostream&
STAHRAeroData::Restart(std::ostream& out) const
{
	switch (profile) {
	case 1:
		out << "NACA0012";
		break;

	case 2:
		out << "RAE9671";
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return RestartUnsteady(out);
}

int
STAHRAeroData::GetForces(int i, const doublereal* W, doublereal* TNG,
		outa_t& OUTA)
{
	switch (unsteadyflag) {
	case AeroData::HARRIS:
	case AeroData::BIELAWA:
		Predict(i, atan2(-W[1], W[0]), OUTA.alf1, OUTA.alf2);
		break;

	default:
		break;
	}

	integer u = unsteadyflag;
	__FC_DECL__(aerod2)(const_cast<doublereal *>(W), reinterpret_cast<doublereal *>(&VAM),
		TNG, reinterpret_cast<doublereal *>(&OUTA), &u, &Omega, &profile);

	return 0;
}

int
STAHRAeroData::GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
	return AeroData::GetForcesJacForwardDiff_int(i, W, TNG, J, OUTA);
}

/* STAHRAeroData - end */

/* C81AeroData - begin */

C81AeroData::C81AeroData(AeroData::UnsteadyModel u, integer p,
		const c81_data* d, DriveCaller *ptime)
: AeroData(u, ptime), profile(p), data(d)
{
	ASSERT(data != NULL);
}

C81AeroData::~C81AeroData(void)
{
	NO_OP;
}

std::ostream&
C81AeroData::Restart(std::ostream& out) const
{
	out << "C81, " << profile;

	return RestartUnsteady(out);
}

int
C81AeroData::GetForces(int i, const doublereal* W, doublereal* TNG, outa_t& OUTA)
{
	switch (unsteadyflag) {
	case AeroData::HARRIS:
	case AeroData::BIELAWA:
		Predict(i, atan2(-W[1], W[0]), OUTA.alf1, OUTA.alf2);
		break;

	default:
		break;
	}

	return c81_aerod2_u(const_cast<doublereal *>(W), &VAM, TNG, &OUTA,
		const_cast<c81_data *>(data), unsteadyflag);
}

int
C81AeroData::GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
	return AeroData::GetForcesJacForwardDiff_int(i, W, TNG, J, OUTA);
}

/* C81AeroData - end */

/* C81MultipleAeroData - begin */

C81MultipleAeroData::C81MultipleAeroData(
		AeroData::UnsteadyModel u,
		integer np,
		integer *p,
		doublereal *ub,
		const c81_data** d,
		DriveCaller *ptime
)
: AeroData(u, ptime), nprofiles(np), profiles(p), upper_bounds(ub), data(d)
{
	ASSERT(nprofiles > 0);
	ASSERT(profiles != NULL);
	ASSERT(upper_bounds != NULL);
	ASSERT(data != NULL);
}

C81MultipleAeroData::~C81MultipleAeroData(void)
{
	SAFEDELETEARR(profiles);
	SAFEDELETEARR(upper_bounds);
	SAFEDELETEARR(data);
}

std::ostream&
C81MultipleAeroData::Restart(std::ostream& out) const
{
	out << "C81, multiple";
	for (int i = 0; i < nprofiles; i++) {
		out << ", " << profiles[i] << ", " << upper_bounds[i];
	}

	return RestartUnsteady(out);
}

void
C81MultipleAeroData::SetSectionData(
		const doublereal& abscissa,
		const doublereal& chord,
		const doublereal& forcepoint,
		const doublereal& velocitypoint,
		const doublereal& twist,
		const doublereal& omega
)
{
	ASSERT(abscissa >= -1. && abscissa <= 1.);

	AeroData::SetSectionData(abscissa, chord, forcepoint, velocitypoint,
			twist, omega);

	for (int i = nprofiles - 1; i--; ) {
		if (abscissa > upper_bounds[i]) {
			curr_data = i + 1;
			return;
		}
	}

	curr_data = 0;
}

int
C81MultipleAeroData::GetForces(int i, const doublereal* W, doublereal* TNG, outa_t& OUTA)
{
	switch (unsteadyflag) {
	case AeroData::HARRIS:
	case AeroData::BIELAWA:
		Predict(i, atan2(-W[1], W[0]), OUTA.alf1, OUTA.alf2);
		break;

	default:
		break;
	}

	return c81_aerod2_u(const_cast<doublereal *>(W), &VAM, TNG, &OUTA,
		const_cast<c81_data *>(data[curr_data]), unsteadyflag);
}

int
C81MultipleAeroData::GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
	return AeroData::GetForcesJacForwardDiff_int(i, W, TNG, J, OUTA);
}

/* C81MultipleAeroData - end */

#ifdef USE_C81INTERPOLATEDAERODATA

/* C81InterpolatedAeroData - begin */

C81InterpolatedAeroData::C81InterpolatedAeroData(
		AeroData::UnsteadyModel u,
		integer np,
		integer *p,
		doublereal *ub,
		const c81_data** d,
		integer i_p,
		DriveCaller *ptime
)
: AeroData(u, ptime), nprofiles(np), profiles(p), upper_bounds(ub), data(d),
i_points(i_p), i_data(0)
{
	ASSERT(nprofiles > 0);
	ASSERT(profiles != NULL);
	ASSERT(upper_bounds != NULL);
	ASSERT(data != NULL);
	ASSERT(i_points > 0);

	SAFENEWARRNOFILL(i_data, c81_data, i_points);

	GaussDataIterator GDI(i_points);
	PntWght PW = GDI.GetFirst();
	do {
		doublereal dCsi = PW.dGetPnt();
		int	from = -1, to;

		for (int i = nprofiles - 1; i--; ) {
			if (dCsi >= upper_bounds[i]) {
				from = i;
				break;
			}
		}

		if (from == -1) {
			silent_cerr("cannot find C81 data lower bound for point xi=" << dCsi << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		ASSERT(from < nprofiles);

		to = from + 1;

		/* we need to interpolate between data[from]
		 * and data[to] */

		/* we only accept homogeneous data sources,
		 * i.e. same Mach and alpha patterns */
		if (data[from]->NML != data[to]->NML) {
			silent_cerr("number of Mach points for Cl between profiles "
					<< from << " (" << data[from]->NML << ") and "
					<< to << " (" << data[to]->NML << ") do not match"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (data[from]->NAL != data[to]->NAL) {
			silent_cerr("number of AoA points for Cl between profiles "
					<< from << " (" << data[from]->NAL << ") and "
					<< to << " (" << data[to]->NAL << ") do not match"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (data[from]->NMD != data[to]->NMD) {
			silent_cerr("number of Mach points for Cd between profiles "
					<< from << " (" << data[from]->NMD << ") and "
					<< to << " (" << data[to]->NMD << ") do not match"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (data[from]->NAD != data[to]->NAD) {
			silent_cerr("number of AoA points for Cd between profiles "
					<< from << " (" << data[from]->NAD << ") and "
					<< to << " (" << data[to]->NAD << ") do not match"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (data[from]->NMM != data[to]->NMM) {
			silent_cerr("number of Mach points for Cm between profiles "
					<< from << " (" << data[from]->NMM << ") and "
					<< to << " (" << data[to]->NMM << ") do not match"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (data[from]->NAM != data[to]->NAM) {
			silent_cerr("number of AoA points for Cm between profiles "
					<< from << " (" << data[from]->NAM << ") and "
					<< to << " (" << data[to]->NAM << ") do not match"
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		for (int i = 0; i < data[from]->NML; i++) {
			if (data[from]->ml[i] != data[to]->ml[i]) {
				silent_cerr("Mach point " << i << "for profiles "
						<< from << " (" << data[from]->ml[i] << ") and "
						<< to << " (" << data[to]->ml[i] << ") differs"
						<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		for (int i = 0; i < data[from]->NAL; i++) {
			if (data[from]->al[i] != data[to]->al[i]) {
				silent_cerr("AoA point " << i << "for profiles "
						<< from << " (" << data[from]->al[i] << ") and "
						<< to << " (" << data[to]->al[i] << ") differs"
						<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

	} while (GDI.fGetNext(PW));
}

C81InterpolatedAeroData::~C81InterpolatedAeroData(void)
{
	SAFEDELETEARR(profiles);
	SAFEDELETEARR(upper_bounds);
	SAFEDELETEARR(data);
}

std::ostream&
C81InterpolatedAeroData::Restart(std::ostream& out) const
{
	out << "C81, interpolated";
	for (int i = 0; i < nprofiles; i++) {
		out << ", " << profiles[i] << ", " << upper_bounds[i];
	}

	return RestartUnsteady(out);
}

void
C81InterpolatedAeroData::SetSectionData(
		const doublereal& abscissa,
		const doublereal& chord,
		const doublereal& forcepoint,
		const doublereal& velocitypoint,
		const doublereal& twist,
		const doublereal& omega
)
{
	ASSERT(abscissa >= -1. && abscissa <= 1.);

	AeroData::SetSectionData(abscissa, chord, forcepoint, velocitypoint,
		twist, omega);

	for (int i = nprofiles - 1; i--; ) {
		if (abscissa > upper_bounds[i]) {
			curr_data = i + 1;
			return;
		}
	}

	curr_data = 0;
}

int
C81InterpolatedAeroData::GetForces(int i, const doublereal* W, doublereal* TNG, outa_t& OUTA)
{
	switch (unsteadyflag) {
	case AeroData::HARRIS:
	case AeroData::BIELAWA: {
		Predict(i, atan2(-W[1], W[0]), OUTA.alf1, OUTA.alf2);
		break;
	}

	default:
		break;
	}

	return c81_aerod2_u(const_cast<doublereal *>(W), &VAM, TNG, &OUTA,
		const_cast<c81_data *>(data[curr_data]), unsteadyflag);
}

int
C81InterpolateAeroData::GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
	return AeroData::GetForcesJacForwardDiff_int(i, W, TNG, J, OUTA);
}

/* C81InterpolatedAeroData - end */

#endif /* USE_C81INTERPOLATEDAERODATA */

/* UMDAeroData - begin */

UMDAeroData::UMDAeroData(DriveCaller *ptime)
: AeroData(STEADY, ptime)
{
	NO_OP;
}

std::ostream&
UMDAeroData::Restart(std::ostream& out) const
{
	return out;
}

// aerodynamic models with internal states
unsigned int
UMDAeroData::iGetNumDof(void) const
{
	return 2;
}

DofOrder::Order
UMDAeroData::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

void
UMDAeroData::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr,
	integer iFirstIndex, integer iFirstSubIndex,
	int i, const doublereal* W, doublereal* TNG, outa_t& OUTA)
{
	// doublereal q1 = XCurr(iFirstIndex + 1);
	// doublereal q2 = XCurr(iFirstIndex + 2);
	doublereal q1p = XPrimeCurr(iFirstIndex + 1);
	doublereal q2p = XPrimeCurr(iFirstIndex + 2);

	doublereal alpha = -atan2(W[1], W[0]);
	doublereal omega = W[5];

	WorkVec.PutCoef(iFirstSubIndex + 1, alpha - q1p);
	WorkVec.PutCoef(iFirstSubIndex + 2, omega - q2p);

	TNG[0] = 0.;
	TNG[1] = 2*M_PI*alpha;
	TNG[2] = 0.;
	TNG[3] = 0.;
	TNG[4] = 0.;
	TNG[5] = 0.;
}

void
UMDAeroData::AssJac(FullSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr,
	integer iFirstIndex, integer iFirstSubIndex,
	const Mat3xN& vx, const Mat3xN& wx, Mat3xN& fq, Mat3xN& cq,
	int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
	// doublereal q1 = XCurr(iFirstIndex + 1);
	// doublereal q2 = XCurr(iFirstIndex + 2);
	// doublereal q1p = XPrimeCurr(iFirstIndex + 1);
	// doublereal q2p = XPrimeCurr(iFirstIndex + 2);

	doublereal alpha = -atan2(W[1], W[0]);
	// doublereal omega = W[5];

	WorkMat.IncCoef(iFirstSubIndex + 1, iFirstSubIndex + 1, 1.);
	WorkMat.IncCoef(iFirstSubIndex + 2, iFirstSubIndex + 2, 1.);

	doublereal dd = W[0]*W[0] + W[1]*W[1];
	doublereal dalpha_dvx = W[1]/dd;
	doublereal dalpha_dvy = -W[0]/dd;
	for (integer iCol = 1; iCol <= 6; iCol++) {
		WorkMat.DecCoef(iFirstSubIndex + 1, iCol,
			dalpha_dvx*vx(1, iCol) + dalpha_dvy*vx(2, iCol));
		WorkMat.DecCoef(iFirstSubIndex + 2, iCol, wx(3, iCol));
	}

	TNG[0] = 0.;
	TNG[1] = 2*M_PI*alpha;
	TNG[2] = 0.;
	TNG[3] = 0.;
	TNG[4] = 0.;
	TNG[5] = 0.;

	J.Reset();
	J(2, 1) = 2*M_PI*dalpha_dvx;
	J(2, 2) = 2*M_PI*dalpha_dvy;
}

/* UMDAeroData - end */

/* C81UnsteadyAeroData - begin */

C81UnsteadyAeroData::C81UnsteadyAeroData(AeroData::UnsteadyModel u, integer p,
		const c81_data* d, DriveCaller *ptime)
: AeroData(u, ptime), profile(p), data(d)
{
	ASSERT(data != NULL);
}

C81UnsteadyAeroData::~C81UnsteadyAeroData(void)
{
	NO_OP;
}

std::ostream&
C81UnsteadyAeroData::Restart(std::ostream& out) const
{
	out << "C81, " << profile;

	return RestartUnsteady(out);
}

// aerodynamic models with internal states
unsigned int
C81UnsteadyAeroData::iGetNumDof(void) const
{
	return 2;
}

DofOrder::Order
C81UnsteadyAeroData::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

void
C81UnsteadyAeroData::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr,
	integer iFirstIndex, integer iFirstSubIndex,
	int i, const doublereal* W, doublereal* TNG, outa_t& OUTA)
{
	doublereal q1p = XPrimeCurr(iFirstIndex + 1);
	doublereal q2p = XPrimeCurr(iFirstIndex + 2);

	// trivial: 1. * delta_dot{q} = -dot{q}

	WorkVec.PutCoef(iFirstSubIndex + 1, -q1p);
	WorkVec.PutCoef(iFirstSubIndex + 2, -q2p);

	c81_aerod2_u(const_cast<doublereal *>(W), &VAM, TNG, &OUTA,
		const_cast<c81_data *>(data), unsteadyflag);
}

void
C81UnsteadyAeroData::AssJac(FullSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr,
	integer iFirstIndex, integer iFirstSubIndex,
	const Mat3xN& vx, const Mat3xN& wx, Mat3xN& fq, Mat3xN& cq,
	int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
	WorkMat.IncCoef(iFirstSubIndex + 1, iFirstSubIndex + 1, 1.);
	WorkMat.IncCoef(iFirstSubIndex + 2, iFirstSubIndex + 2, 1.);

	AeroData::GetForcesJacForwardDiff_int(i, W, TNG, J, OUTA);

	// probably, we need to reset fq, cq
}

/* C81UnsteadyAeroData - end */
