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

/* C81TheodorsenAeroData - begin */

static const doublereal TheodorsenParams[2][4] = {
	{ 0.165, 0.335, 0.0455, 0.3 },
	{ 0.165, 0.335, 0.041, 0.32 }
};

C81TheodorsenAeroData::C81TheodorsenAeroData(integer p,
	const c81_data* d,
	doublereal omegaPD)
: C81AeroData(STEADY, p, d, 0), iParam(0), omegaPD(omegaPD)
{
	NO_OP;
}

C81TheodorsenAeroData::~C81TheodorsenAeroData(void)
{
	NO_OP;
}

std::ostream&
C81TheodorsenAeroData::Restart(std::ostream& out) const
{
	out << "theodorsen, " << profile;

	return RestartUnsteady(out);
}

// aerodynamic models with internal states
unsigned int
C81TheodorsenAeroData::iGetNumDof(void) const
{
	return 6;
}

DofOrder::Order
C81TheodorsenAeroData::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

void
C81TheodorsenAeroData::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr,
	integer iFirstIndex, integer iFirstSubIndex,
	int i, const doublereal* W, doublereal* TNG, outa_t& OUTA)
{
	doublereal q1 = XCurr(iFirstIndex + 1);
	doublereal q2 = XCurr(iFirstIndex + 2);
	doublereal q3 = XCurr(iFirstIndex + 3);
	doublereal q4 = XCurr(iFirstIndex + 4);
	doublereal q5 = XCurr(iFirstIndex + 5);
	doublereal q6 = XCurr(iFirstIndex + 6);

	doublereal q1p = XPrimeCurr(iFirstIndex + 1);
	doublereal q2p = XPrimeCurr(iFirstIndex + 2);
	doublereal q3p = XPrimeCurr(iFirstIndex + 3);
	doublereal q4p = XPrimeCurr(iFirstIndex + 4);
	doublereal q5p = XPrimeCurr(iFirstIndex + 5);
	doublereal q6p = XPrimeCurr(iFirstIndex + 6);

	doublereal d14 = VAM.force_position;
	doublereal d34 = VAM.bc_position;
	doublereal chord = VAM.chord;

	doublereal a = (d14 + d34)/chord;

	doublereal Uinf = sqrt(W[0]*W[0] + W[1]*W[1]);

	doublereal A1 = TheodorsenParams[iParam][0];
	doublereal A2 = TheodorsenParams[iParam][1];
	doublereal b1 = TheodorsenParams[iParam][2];
	doublereal b2 = TheodorsenParams[iParam][3];

	doublereal u1 = atan2(- W[1] - W[5]*d14, W[0]);
	doublereal u2 = atan2(- W[1] - W[5]*d34, W[0]);

	doublereal d = 2*Uinf/chord;

	doublereal y1 = (A1 + A2)*b1*b2*d*d*q1 + (A1*b1 + A2*b2)*d*q2
		+ (1 - A1 - A2)*u2;
	doublereal y2 = omegaPD*omegaPD*q3;
	doublereal y3 = omegaPD*omegaPD*q5;
	doublereal y4 = u2/((d34 - d14)*Uinf);

	doublereal tan_y1 = std::tan(y1);
	doublereal Vxp2 = Uinf*Uinf - pow(W[0]*tan_y1, 2)
		- std::pow(d34*W[5], 2);

	doublereal WW[6];
	WW[0] = copysign(std::sqrt(Vxp2), W[0]);
	WW[1] = -W[0]*tan_y1 - d34*W[5];
	WW[2] = W[2];
	WW[3] = W[3];
	WW[4] = W[4];
	WW[5] = W[5];

	c81_aerod2_u(WW, &VAM, TNG, &OUTA,
		const_cast<c81_data *>(data), unsteadyflag);

	doublereal UUinf2 = Uinf*Uinf + W[2]*W[2];
	doublereal qD = .5*VAM.density*UUinf2;

	doublereal clalpha = OUTA.clalpha;
	if (clalpha > 0.) {
		TNG[1] += qD*chord*clalpha/2/d*(y2 - a/d*y3);
		TNG[5] += qD*chord*chord*clalpha/2/d*(-y2/4 + (a*a/d + a/4/d - 1/d/16)*y3 - y4/4);
	}

	WorkVec.PutCoef(iFirstSubIndex + 1, -q1p + q2);
	WorkVec.PutCoef(iFirstSubIndex + 2, -q2p - b1*b2*d*d*q1 - (b1 + b2)*d*q2 + u2);
	WorkVec.PutCoef(iFirstSubIndex + 3, -q3p - 2*omegaPD*q3 - omegaPD*omegaPD*q4 + (d34*u1 - d14*u2)/(d34 - d14));
	WorkVec.PutCoef(iFirstSubIndex + 4, -q4p + q3);
	WorkVec.PutCoef(iFirstSubIndex + 5, -q5p - 2*omegaPD*q5 - omegaPD*omegaPD*q6 - (u1 - u2)/(d34 - d14)/Uinf);
	WorkVec.PutCoef(iFirstSubIndex + 6, -q6p + q5);
}

void
C81TheodorsenAeroData::AssJac(FullSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr,
	integer iFirstIndex, integer iFirstSubIndex,
	const Mat3xN& vx, const Mat3xN& wx, Mat3xN& fq, Mat3xN& cq,
	int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
	// doublereal A1 = TheodorsenParams[iParam][0];
	// doublereal A2 = TheodorsenParams[iParam][1];
	doublereal b1 = TheodorsenParams[iParam][2];
	doublereal b2 = TheodorsenParams[iParam][3];

	doublereal Uinf = sqrt(W[0]*W[0] + W[1]*W[1]);
	doublereal d = 2*Uinf/VAM.chord;

	WorkMat.IncCoef(iFirstSubIndex + 1, iFirstSubIndex + 1, 1.);
	WorkMat.IncCoef(iFirstSubIndex + 2, iFirstSubIndex + 2, 1.);
	WorkMat.IncCoef(iFirstSubIndex + 3, iFirstSubIndex + 3, 1.);
	WorkMat.IncCoef(iFirstSubIndex + 4, iFirstSubIndex + 4, 1.);
	WorkMat.IncCoef(iFirstSubIndex + 5, iFirstSubIndex + 5, 1.);
	WorkMat.IncCoef(iFirstSubIndex + 6, iFirstSubIndex + 6, 1.);

	WorkMat.IncCoef(iFirstSubIndex + 1, iFirstSubIndex + 2, -dCoef);
	WorkMat.IncCoef(iFirstSubIndex + 2, iFirstSubIndex + 1, dCoef*b1*b2*d*d);
	WorkMat.IncCoef(iFirstSubIndex + 2, iFirstSubIndex + 2, dCoef*(b1 + b2)*d);

	WorkMat.IncCoef(iFirstSubIndex + 3, iFirstSubIndex + 3, dCoef*2*omegaPD);
	WorkMat.IncCoef(iFirstSubIndex + 3, iFirstSubIndex + 4, dCoef*omegaPD*omegaPD);
	WorkMat.IncCoef(iFirstSubIndex + 4, iFirstSubIndex + 3, -dCoef);

	WorkMat.IncCoef(iFirstSubIndex + 5, iFirstSubIndex + 5, dCoef*2*omegaPD);
	WorkMat.IncCoef(iFirstSubIndex + 5, iFirstSubIndex + 6, dCoef*omegaPD*omegaPD);
	WorkMat.IncCoef(iFirstSubIndex + 6, iFirstSubIndex + 5, -dCoef);

	AeroData::GetForcesJacForwardDiff_int(i, W, TNG, J, OUTA);

	// probably, we need to reset fq, cq
}

/* C81TheodorsenAeroData - end */
