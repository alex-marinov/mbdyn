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
		Predict(i, atan2(-W[VY], W[VX]), OUTA.alf1, OUTA.alf2);
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
		Predict(i, atan2(-W[VY], W[VX]), OUTA.alf1, OUTA.alf2);
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
		Predict(i, atan2(-W[VY], W[VX]), OUTA.alf1, OUTA.alf2);
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
		Predict(i, atan2(-W[VY], W[VX]), OUTA.alf1, OUTA.alf2);
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

static const doublereal TheodorsenParams[2][4] = {
	{ 0.165, 0.335, 0.0455, 0.3 },
	{ 0.165, 0.335, 0.041, 0.32 }
};

#ifdef USE_UNSTEADY_6STATES

/* C81TheodorsenAeroData - begin */


C81TheodorsenAeroData::C81TheodorsenAeroData(integer p,
	const c81_data* d,
	doublereal omegaPD, DriveCaller *ptime)
: C81AeroData(STEADY, p, d, ptime), iParam(0), omegaPD(omegaPD)
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

	d14 = VAM.force_position;
	d34 = VAM.bc_position;
	chord = VAM.chord;

	a = (d14 + d34)/chord;

	doublereal Uinf = sqrt(W[VX]*W[VX] + W[VY]*W[VY]);

	A1 = TheodorsenParams[iParam][0];
	A2 = TheodorsenParams[iParam][1];
	b1 = TheodorsenParams[iParam][2];
	b2 = TheodorsenParams[iParam][3];

	doublereal u1 = atan2(- W[VY] - W[WZ]*d14, W[VX]);
	doublereal u2 = atan2(- W[VY] - W[WZ]*d34, W[VX]);

	doublereal d = 2*Uinf/chord;

	doublereal y1 = (A1 + A2)*b1*b2*d*d*q1 + (A1*b1 + A2*b2)*d*q2
		+ (1 - A1 - A2)*u2;
	doublereal y2 = omegaPD*omegaPD*q3;
	doublereal y3 = omegaPD*omegaPD*q5;
	doublereal y4 = (u1-u2)/((d34 - d14)/Uinf);

	doublereal tan_y1 = std::tan(y1);
	doublereal Vxp2 = Uinf*Uinf - pow(W[VX]*tan_y1, 2)
		- std::pow(d34*W[WZ], 2);

	doublereal WW[6];
	WW[VX] = copysign(std::sqrt(Vxp2), W[VX]);
	WW[VY] = -W[VX]*tan_y1 - d34*W[WZ];
	WW[VZ] = W[VZ];
	WW[WX] = W[WX];
	WW[WY] = W[WY];
	WW[WZ] = W[WZ];

	c81_aerod2_u(WW, &VAM, TNG, &OUTA,
		const_cast<c81_data *>(data), unsteadyflag);

	doublereal UUinf2 = Uinf*Uinf + W[VZ]*W[VZ];
	doublereal qD = .5*VAM.density*UUinf2;
	if (qD > std::numeric_limits<doublereal>::epsilon()) {
		cx_0 = TNG[0]/(qD*chord);
		cy_0 = TNG[1]/(qD*chord);
		cz_0 = TNG[2]/(qD*chord);
		cmx_0 = TNG[3]/(qD*chord*chord);
		cmy_0 = TNG[4]/(qD*chord*chord);
		cmz_0 = TNG[5]/(qD*chord*chord);
	}else{
		cx_0 = 0.;
		cy_0 = 0.;
		cz_0 = 0.;
		cmx_0 = 0.;
		cmy_0 = 0.;
		cmz_0 = 0.;
	}


	clalpha = OUTA.clalpha;
	if (clalpha > 0.) {
		TNG[1] += qD*chord*clalpha/2/d*(y2 - a/d*y3);
		TNG[5] += qD*chord*chord*clalpha/2/d*(-y2/4 + (a/4/d - 1/d/16)*y3 - y4/4);
	}

	WorkVec.PutCoef(iFirstSubIndex + 1, -q1p + q2);
	WorkVec.PutCoef(iFirstSubIndex + 2, -q2p - b1*b2*d*d*q1 - (b1 + b2)*d*q2 + u2);
	WorkVec.PutCoef(iFirstSubIndex + 3, -q3p - 2*omegaPD*q3 - omegaPD*omegaPD*q4 + (d34*u1 - d14*u2)/(d34 - d14));
	WorkVec.PutCoef(iFirstSubIndex + 4, -q4p + q3);
	/* mbdyn e theodorsen hanno l'asse x in direzione opposta e
 	* e quindi d14 e d34 sono di segno opposto! */
	WorkVec.PutCoef(iFirstSubIndex + 5, -q5p - 2*omegaPD*q5 - omegaPD*omegaPD*q6 + (u1 - u2)/(d34 - d14)*Uinf);
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

	doublereal q1 = XCurr(iFirstIndex + 1);
	doublereal q2 = XCurr(iFirstIndex + 2);
	doublereal q3 = XCurr(iFirstIndex + 3);
	doublereal q4 = XCurr(iFirstIndex + 4);
	doublereal q5 = XCurr(iFirstIndex + 5);
	doublereal q6 = XCurr(iFirstIndex + 6);

	doublereal Uinf = sqrt(W[VX]*W[VX] + W[VY]*W[VY]);
	doublereal d = 2*Uinf/chord;
	doublereal UUinf2 = Uinf*Uinf + W[VZ]*W[VZ];
	doublereal qD = .5*VAM.density*UUinf2;

	doublereal u1 = atan2(- W[VY] - W[WZ]*d14, W[VX]);
	doublereal u2 = atan2(- W[VY] - W[WZ]*d34, W[VX]);

	doublereal y1 = (A1 + A2)*b1*b2*d*d*q1 + (A1*b1 + A2*b2)*d*q2
		+ (1 - A1 - A2)*u2;
	doublereal y2 = omegaPD*omegaPD*q3;
	doublereal y3 = omegaPD*omegaPD*q5;
	doublereal y4 = (u1-u2)/((d34 - d14)/Uinf);


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



	/* computing the matrix g_{/\tilde{v}} - dimension: 6x3, where 6 is n_states */
	doublereal dU_V_11, dU_V_12, dU_V_21, dU_V_22;

	doublereal dDen14 = W[VX]*W[VX] + W[VY]*W[VY] + W[WZ]*W[WZ]*d14*d14 + 2*W[VY]*W[WZ]*d14;
	doublereal dDen34 = W[VX]*W[VX] + W[VY]*W[VY] + W[WZ]*W[WZ]*d34*d34 + 2*W[VY]*W[WZ]*d34;

	dU_V_11 = (W[VY] + W[WZ]*d14)/dDen14;
	dU_V_12 = -W[VX]/dDen14;
	dU_V_21 = (W[VY] + W[WZ]*d34)/dDen34;
	dU_V_22 = -W[VX]/dDen34;

	doublereal B31, B32, B51, B52;

	B31 = d34/(d34-d14);
	B32 = -d14/(d34-d14);
	B51 = Uinf/(d34-d14);
	B52 = -Uinf/(d34-d14);


	doublereal dG_V_21, dG_V_22, dG_V_31, dG_V_32, dG_V_51, dG_V_52;

	dG_V_21 =  (b1*b2*d*4*q1/chord + (b1+b2)*2*q2/chord)*W[VX]/Uinf - dU_V_21;
	dG_V_22 =  (b1*b2*d*4*q1/chord + (b1+b2)*2*q2/chord)*W[VY]/Uinf - dU_V_22;
	dG_V_31 = -B31*dU_V_11 - B32*dU_V_21;
	dG_V_32 = -B31*dU_V_12 - B32*dU_V_22;
	dG_V_51 = -( (u1-u2)/(d34-d14) )*W[VX]/Uinf - B51*dU_V_11 - B52*dU_V_21;
	dG_V_52 = -( (u1-u2)/(d34-d14) )*W[VY]/Uinf - B51*dU_V_12 - B52*dU_V_22;

	
	/* computing the matrix g_{/\tilde{\omega}} - dimension: 6x3, where 6 is n_states */
	doublereal dU_W_13, dU_W_23;

	dU_W_13 = - W[VX]*d14/dDen14;
	dU_W_23 = - W[VX]*d34/dDen34;

	doublereal dG_W_23, dG_W_33, dG_W_53;

	dG_W_23 = -dU_W_23;
	dG_W_33 = -B31*dU_W_13 - B32*dU_W_23;
	dG_W_53 = -B51*dU_W_13 - B52*dU_W_23;

	/* assembling the jacobian */
	unsigned int iIndexColumn;
	/* faccio i calcoli tenendo conto della sparsità della matrice g_{\/\tilde{v}} */
	for( iIndexColumn = 1; iIndexColumn <= vx.iGetNumCols(); iIndexColumn++){
		WorkMat.IncCoef(iFirstSubIndex+2, iIndexColumn, dG_V_21*vx(1,iIndexColumn) + dG_V_22*vx(2,iIndexColumn));
		WorkMat.IncCoef(iFirstSubIndex+2, iIndexColumn, dG_W_23*wx(3,iIndexColumn));
		WorkMat.IncCoef(iFirstSubIndex+3, iIndexColumn, dG_V_31*vx(1,iIndexColumn) + dG_V_32*vx(2,iIndexColumn));
		WorkMat.IncCoef(iFirstSubIndex+3, iIndexColumn, dG_W_33*wx(3,iIndexColumn));
		WorkMat.IncCoef(iFirstSubIndex+5, iIndexColumn, dG_V_51*vx(1,iIndexColumn) + dG_V_52*vx(2,iIndexColumn));
		WorkMat.IncCoef(iFirstSubIndex+5, iIndexColumn, dG_W_53*wx(3,iIndexColumn));
	}

	/* computing the matrix fq (f_a_{/q}) 3x6 */
	
	doublereal C11, C12, C23, C35;

	C11 = (A1+A2)*b1*b2*d*d;
	C12 = (A1*b1+A2*b2)*d;
	C23 = omegaPD*omegaPD;
	C35 = omegaPD*omegaPD;

	/* computing the derivative of the aerodynamic coefficient in the lookup table
 	 * using the finite difference method */
	/* perturbazione per calcolo delle derivare con le differenze finite */
	doublereal dDeltaY1 = 0.1*M_PI/180.;
	doublereal y1d = y1 + dDeltaY1;
	doublereal tan_y1 = std::tan(y1d);
	doublereal Vxp2 = Uinf*Uinf - pow(W[VX]*tan_y1, 2) - std::pow(d34*W[WZ], 2);

	doublereal WW[6];
	WW[VX] = copysign(std::sqrt(Vxp2), W[VX]);
	WW[VY] = -W[VX]*tan_y1 - d34*W[WZ];
	WW[VZ] = W[VZ];
	WW[WX] = W[WX];
	WW[WY] = W[WY];
	WW[WZ] = W[WZ];

	c81_aerod2_u(WW, &VAM, TNG, &OUTA, const_cast<c81_data *>(data), unsteadyflag);
	doublereal cx_1, cy_1, cmz_1;
	if (qD > std::numeric_limits<doublereal>::epsilon()) {
		cx_1 = TNG[0]/(qD*chord);
		cy_1 = TNG[1]/(qD*chord);
		cmz_1 = TNG[5]/(qD*chord*chord);
	}else{
		cx_1 = 0.;
		cy_1 = 0.;
		cmz_1 = 0.;
	}

	doublereal dCd_alpha = (cx_1-cx_0)/dDeltaY1;
	doublereal dCl_alpha = (cy_1-cy_0)/dDeltaY1;
	doublereal dCm_alpha = (cmz_1-cmz_0)/dDeltaY1;

	doublereal qDc = qD*chord;

	fq.Put(1, 1, qDc*dCd_alpha*C11);
	fq.Put(1, 2, qDc*dCd_alpha*C12);
	fq.Put(2, 1, qDc*dCl_alpha*C11);
	fq.Put(2, 2, qDc*dCl_alpha*C12);
	fq.Put(2, 3, qDc*0.5*clalpha*C23/d);
	fq.Put(2, 5, -qDc*0.5*clalpha*a*C35/(d*d));
	
	cq.Put(3, 1, qDc*chord*dCm_alpha*C11);
	cq.Put(3, 2, qDc*chord*dCm_alpha*C12);
	cq.Put(3, 3, -qDc*chord*0.25*0.5*clalpha*C23/d);
	cq.Put(3, 5, qDc*chord*0.5*clalpha*( (0.25*a/d) - (1/(16*d)) )*C35/d);

	/* computing the J matrix */
	/* f_{/\tilde{v}} */
	doublereal dY_V_11, dY_V_12, dY_V_41, dY_V_42;
	doublereal D12, D41, D42;

	D12 = 1-A1-A2;
	D41 = Uinf/(d34-d14);
	D42 = -Uinf/(d34-d14);

	/* (C_{/Uinf}*q + D_{/Uinf}*u)*Uinf_{v} */
	dY_V_11 = ((A1+A2)*b1*b2*d*(4./chord)*q1 + (A1*b1+A2*b2)*(2./chord)*q2)*W[VX]/Uinf;
	dY_V_12 = ((A1+A2)*b1*b2*d*(4./chord)*q1 + (A1*b1+A2*b2)*(2./chord)*q2)*W[VY]/Uinf;
	dY_V_41 = ((u1-u2)/(d34-d14))*W[VX]/Uinf;
	dY_V_42 = ((u1-u2)/(d34-d14))*W[VY]/Uinf;
	/* (C_{/Uinf}*q + D_{/Uinf}*u)*Uinf_{v} + D*u_{/v} */
	dY_V_11 += D12*dU_V_21;
	dY_V_12 += D12*dU_V_22;
	dY_V_41 += D41*dU_V_11 + D42*dU_V_21;
	dY_V_42 += D41*dU_V_12 + D42*dU_V_22;

	doublereal dCfy_Uinf = 0.5*clalpha*(-y2 + (chord*a*y3)/Uinf)/(d*Uinf);
	doublereal rho = VAM.density;

	doublereal cy = cy_0 + clalpha/2/d*(y2 - a/d*y3);

	J.Put(1, 1, rho*chord*cx_0*W[VX] + qDc*( dCd_alpha*dY_V_11));
	J.Put(1, 2, rho*chord*cx_0*W[VY] + qDc*( dCd_alpha*dY_V_12));
	J.Put(1, 3, rho*chord*cx_0*W[VZ]);
	J.Put(2, 1, rho*chord*cy*W[VX] + qDc*( dCl_alpha*dY_V_11 + dCfy_Uinf*W[VX]/Uinf) );
	J.Put(2, 2, rho*chord*cy*W[VY] + qDc*( dCl_alpha*dY_V_12 + dCfy_Uinf*W[VY]/Uinf ) );
	J.Put(2, 3, rho*chord*cy*W[VZ]);
	J.Put(3, 1, rho*chord*cz_0*W[VX]);
	J.Put(3, 2, rho*chord*cz_0*W[VY]);
	J.Put(3, 3, rho*chord*cz_0*W[VZ]);

	/* f_{/\tilde{omega}} */
	J.Put(1, 6, qDc*dCd_alpha*D12*dU_W_23);
	J.Put(2, 6, qDc*dCl_alpha*D12*dU_W_23);

	/* c_{/\tilde{v}} */
	doublereal dCmz_Uinf = 0.5*clalpha*(y2 - ( ((chord*a)/(Uinf)) - ((chord)/(4*Uinf)) )*y3 +y4)/(d*4*Uinf);

	doublereal cmz = cmz_0 + clalpha/2/d*( -y2/4 + (a/4/d -1/d/16)*y3 -y4/4);


	J.Put(4, 1, rho*chord*chord*cmx_0*W[VX]);
	J.Put(4, 2, rho*chord*chord*cmx_0*W[VY]);
	J.Put(4, 3, rho*chord*chord*cmx_0*W[VZ]);
	J.Put(5, 1, rho*chord*chord*cmy_0*W[VX]);
	J.Put(5, 2, rho*chord*chord*cmy_0*W[VY]);
	J.Put(5, 3, rho*chord*chord*cmy_0*W[VZ]);
	J.Put(6, 1, rho*chord*chord*cmz*W[VX] +qDc*chord*( dCm_alpha*dY_V_11 + (-0.25*0.5*clalpha/d)*dY_V_41 + dCmz_Uinf*W[VX]/Uinf ) );
	J.Put(6, 2, rho*chord*chord*cmz*W[VY] +qDc*chord*( dCm_alpha*dY_V_12 + (-0.25*0.5*clalpha/d)*dY_V_42 + dCmz_Uinf*W[VY]/Uinf ) );
	J.Put(6, 3, rho*chord*chord*cmz*W[VZ]);

	/* c_{/\tilde{omega}} */
	J.Put(6, 6, qDc*chord*( dCm_alpha*D12*dU_W_23 + (-0.25*0.5*clalpha/d)*(D41*dU_W_13 + D42*dU_W_23) ) );

#ifdef  DEBUG_JACOBIAN_UNSTEADY	

	printf("G/v matrix\n");
	printf("%lf %lf %lf\n", 0., 0., 0.);
	printf("%lf %lf %lf\n", dG_V_21, dG_V_22, 0.);
	printf("%lf %lf %lf\n", dG_V_31, dG_V_32, 0.);
	printf("%lf %lf %lf\n", 0., 0., 0.);
	printf("%lf %lf %lf\n", dG_V_51, dG_V_52, 0.);
	printf("%lf %lf %lf\n", 0., 0., 0.);

	printf("G/w matrix\n");
	printf("%lf %lf %lf\n", 0., 0., 0.);
	printf("%lf %lf %lf\n", 0., 0., dG_W_23);
	printf("%lf %lf %lf\n", 0., 0., dG_W_33);
	printf("%lf %lf %lf\n", 0., 0., 0.);
	printf("%lf %lf %lf\n", 0., 0., dG_W_53);
	printf("%lf %lf %lf\n", 0., 0., 0.);

	printf("f/q matrix\n");
	printf("%lf %lf %lf %lf %lf %lf\n", fq(1,1), fq(1,2), 0., 0., 0., 0.);
	printf("%lf %lf %lf %lf %lf %lf\n", fq(2,1), fq(2,2), fq(2,3), 0., fq(2,5), 0.);
	printf("%lf %lf %lf %lf %lf %lf\n", 0., 0., 0., 0., 0., 0.);

	printf("c/q matrix\n");
	printf("%lf %lf %lf %lf %lf %lf\n", 0., 0., 0., 0., 0., 0.);
	printf("%lf %lf %lf %lf %lf %lf\n", 0., 0., 0., 0., 0., 0.);
	printf("%lf %lf %lf %lf %lf %lf\n", cq(3,1), cq(3,2), cq(3,3), 0., cq(3,5), 0.);

	printf("J matrix\n");
	doublereal iii, jjj;
	for( iii=1; iii<=6; iii++){
		for( jjj=1; jjj<=6; jjj++){
			printf("%lf ", J(iii,jjj));
		}
		printf("\n");
	}
	
	FILE *fd;
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/X.mat","w");
	printf("DATI X\n");
	for( iii=1; iii<=6; iii++){
		printf("%lf ", XCurr(iFirstIndex+iii));
		fprintf(fd,"%15.7e ", XCurr(iFirstIndex+iii));
	}
	fclose(fd);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/XP.mat","w");
	printf("\n");
	printf("DATI XP\n");
	for( iii=1; iii<=6; iii++){
		printf("%lf ", XPrimeCurr(iFirstIndex+iii));
		fprintf(fd,"%15.7e ", XPrimeCurr(iFirstIndex+iii));
	}
	printf("\n");
	fclose(fd);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/W.mat","w");
	printf("DATI W\n");
	printf("%lf %lf %lf %lf %lf %lf", W[VX], W[VY], W[VZ], W[WX], W[WY], W[WZ]);
	fprintf(fd,"%15.7e %15.7e %15.7e %15.7e %15.7e %15.7e", W[VX], W[VY], W[VZ], W[WX], W[WY], W[WZ]);
	fclose(fd);
	printf("\n");
	printf("PARAMETRI \n");
	printf("d14 %lf ", d14);
	printf("\nd34 %lf ", d34);
	printf("\nA1 %lf ", A1);
	printf("\nA2 %lf ", A2);
	printf("\nb1 %lf ", b1);
	printf("\nb2 %lf ", b2);
	printf("\nchord %lf ", chord);
	printf("\na %lf ", a);
	printf("\nrho %lf ", rho);
	printf("\nomegaPD %lf ", omegaPD);
	printf("\ncx_0 %lf ", cx_0);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/cx_0.mat","w");
	fprintf(fd,"%15.7e", cx_0);
	fclose(fd);
	printf("\ncy_0 %lf ", cy_0);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/cy_0.mat","w");
	fprintf(fd,"%15.7e ", cy_0);
	fclose(fd);
	printf("\ncz_0 %lf ", cz_0);
	printf("\ncmx_0 %lf ", cmx_0);
	printf("\ncmy_0 %lf ", cmy_0);
	printf("\ncmz_0 %lf ", cmz_0);
	printf("\ndCl_alpha %lf ", dCl_alpha);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/dCl_alpha.mat","w");
	fprintf(fd,"%15.7e", dCl_alpha);
	fclose(fd);
	printf("\ndCd_alpha %lf", dCd_alpha);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/dCd_alpha.mat","w");
	fprintf(fd,"%15.7e", dCd_alpha);
	fclose(fd);
	printf("\ndCm_alpha %lf", dCm_alpha);
	printf("\nclalpha %lf ", clalpha);
	printf("\n q1 q2 i %lf %lf %d", q1, q2, i);
	getchar();	
#endif /* DEBUG_JACOBIAN_UNSTEADY */	



	AeroData::GetForcesJacForwardDiff_int(i, W, TNG, J, OUTA);

	// probably, we need to reset fq, cq
}

/* C81TheodorsenAeroData - end */

#endif /* USE_UNSTEADY_6STATES */

#ifdef USE_UNSTEADY_2STATES

/* C81TheodorsenAeroData2 - begin */
C81TheodorsenAeroData::C81TheodorsenAeroData(integer p,
	const c81_data* d,
	doublereal omegaPD, DriveCaller *ptime)
: C81AeroData(STEADY, p, d, ptime), iParam(0), omegaPD(omegaPD)
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
	return 2;
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

	doublereal q1p = XPrimeCurr(iFirstIndex + 1);
	doublereal q2p = XPrimeCurr(iFirstIndex + 2);

	d14 = VAM.force_position;
	d34 = VAM.bc_position;
	chord = VAM.chord;

	a = (d14 + d34)/chord;

	doublereal Uinf = sqrt(W[VX]*W[VX] + W[VY]*W[VY]);

	A1 = TheodorsenParams[iParam][0];
	A2 = TheodorsenParams[iParam][1];
	b1 = TheodorsenParams[iParam][2];
	b2 = TheodorsenParams[iParam][3];

	doublereal u1 = atan2(- W[VY] - W[WZ]*d14, W[VX]);
	doublereal u2 = atan2(- W[VY] - W[WZ]*d34, W[VX]); //deve essere l'incidenza a 3/4 corda!!!

	doublereal d = 2*Uinf/chord;

	doublereal y1 = (A1 + A2)*b1*b2*d*d*q1 + (A1*b1 + A2*b2)*d*q2
		+ (1 - A1 - A2)*u2;

	doublereal tan_y1 = std::tan(y1);
	doublereal Vxp2 = Uinf*Uinf - pow(W[VX]*tan_y1, 2)
		- std::pow(d34*W[WZ], 2);

	doublereal WW[6];
	WW[VX] = copysign(std::sqrt(Vxp2), W[VX]);
	WW[VY] = -W[VX]*tan_y1 - d34*W[WZ];
	WW[VZ] = W[VZ];
	WW[WX] = W[WX];
	WW[WY] = W[WY];
	WW[WZ] = W[WZ];

	c81_aerod2_u(WW, &VAM, TNG, &OUTA,
		const_cast<c81_data *>(data), unsteadyflag);

	doublereal UUinf2 = Uinf*Uinf + W[VZ]*W[VZ];
	doublereal qD = .5*VAM.density*UUinf2;
	if (qD > std::numeric_limits<doublereal>::epsilon()) {
		cx_0 = TNG[0]/(qD*chord);
		cy_0 = TNG[1]/(qD*chord);
		cz_0 = TNG[2]/(qD*chord);
		cmx_0 = TNG[3]/(qD*chord*chord);
		cmy_0 = TNG[4]/(qD*chord*chord);
		cmz_0 = TNG[5]/(qD*chord*chord);
	}else{
		cx_0 = 0.;
		cy_0 = 0.;
		cz_0 = 0.;
		cmx_0 = 0.;
		cmy_0 = 0.;
		cmz_0 = 0.;
	}

	alpha_pivot = (u1*d34 - u2*d14)/(d34-d14);
	dot_alpha = Uinf*((u1 - u2)/(d34-d14));

	doublereal Delta_t = pTime->dGet() - prev_time;
	/* controllo che dovrebbe servire solo al primo istante di
 	 * di tempo, in cui non calcolo le derivate per differenze
 	 * finite per non avere 0 a denominatore. */ 
	if (Delta_t > std::numeric_limits<doublereal>::epsilon()) {
		dot_alpha_pivot = (alpha_pivot-prev_alpha_pivot)/Delta_t;
		ddot_alpha = (dot_alpha-prev_dot_alpha)/Delta_t;
	} else {
		printf("Delta time %lf\n", Delta_t);
		dot_alpha_pivot = 0.;
		ddot_alpha = 0.;
	}

	clalpha = OUTA.clalpha;
	if (clalpha > 0.) {
		TNG[1] += qD*chord*clalpha/2/d*(dot_alpha_pivot - a/d*ddot_alpha);
		TNG[5] += qD*chord*chord*clalpha/2/d*(-dot_alpha_pivot/4 + (a/4/d - 1/d/16)*ddot_alpha - dot_alpha/4);
	}

	WorkVec.PutCoef(iFirstSubIndex + 1, -q1p + q2);
	WorkVec.PutCoef(iFirstSubIndex + 2, -q2p - b1*b2*d*d*q1 - (b1 + b2)*d*q2 + u2);
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

	doublereal q1 = XCurr(iFirstIndex + 1);
	doublereal q2 = XCurr(iFirstIndex + 2);

	doublereal Uinf = sqrt(W[VX]*W[VX] + W[VY]*W[VY]);
	doublereal d = 2*Uinf/chord;
	doublereal UUinf2 = Uinf*Uinf + W[VZ]*W[VZ];
	doublereal qD = .5*VAM.density*UUinf2;

	doublereal u1 = atan2(- W[VY] - W[WZ]*d14, W[VX]);
	doublereal u2 = atan2(- W[VY] - W[WZ]*d34, W[VX]);

	doublereal y1 = (A1 + A2)*b1*b2*d*d*q1 + (A1*b1 + A2*b2)*d*q2
		+ (1 - A1 - A2)*u2;


	WorkMat.IncCoef(iFirstSubIndex + 1, iFirstSubIndex + 1, 1.);
	WorkMat.IncCoef(iFirstSubIndex + 2, iFirstSubIndex + 2, 1.);

	WorkMat.IncCoef(iFirstSubIndex + 1, iFirstSubIndex + 2, -dCoef);
	WorkMat.IncCoef(iFirstSubIndex + 2, iFirstSubIndex + 1, dCoef*b1*b2*d*d);
	WorkMat.IncCoef(iFirstSubIndex + 2, iFirstSubIndex + 2, dCoef*(b1 + b2)*d);


	/* computing the matrix g_{/\tilde{v}} - dimension: 2x3, where 2 is n_states */
	doublereal dU_V_11, dU_V_12;

	doublereal dDen34 = W[VX]*W[VX] + W[VY]*W[VY] + W[WZ]*W[WZ]*d34*d34 + 2*W[VY]*W[WZ]*d34;

	dU_V_11 = (W[VY] + W[WZ]*d34)/dDen34;
	dU_V_12 = -W[VX]/dDen34;

	doublereal dG_V_21, dG_V_22;

	dG_V_21 =  (b1*b2*d*4*q1/chord + (b1+b2)*2*q2/chord)*W[VX]/Uinf - dU_V_11;
	dG_V_22 =  (b1*b2*d*4*q1/chord + (b1+b2)*2*q2/chord)*W[VY]/Uinf - dU_V_12;

	
	/* computing the matrix g_{/\tilde{\omega}} - dimension: 2x3, where 2 is n_states */
	doublereal dU_W_13;

	dU_W_13 = - W[VX]*d34/dDen34;

	doublereal dG_W_23;

	dG_W_23 = -dU_W_13;

	/* assembling the jacobian */
	unsigned int iIndexColumn;
	/* faccio i calcoli tenendo conto della sparsità della matrice g_{\/\tilde{v}} */
	for( iIndexColumn = 1; iIndexColumn <= vx.iGetNumCols(); iIndexColumn++){
		WorkMat.IncCoef(iFirstSubIndex+2, iIndexColumn, dG_V_21*vx(1,iIndexColumn) + dG_V_22*vx(2,iIndexColumn));
		WorkMat.IncCoef(iFirstSubIndex+2, iIndexColumn, dG_W_23*wx(3,iIndexColumn));
	}
	/* computing the matrix fq (f_a_{/q}) 3x2 */
	
	/* computing the derivative of the aerodynamic coefficient in the lookup table
 	 * using the finite difference method */
	/* perturbazione per calcolo delle derivare con le differenze finite */
	doublereal dDeltaY1 = 0.1*M_PI/180.;
	doublereal y1d = y1 + dDeltaY1;
	doublereal tan_y1 = std::tan(y1d);
	doublereal Vxp2 = Uinf*Uinf - pow(W[VX]*tan_y1, 2) - std::pow(d34*W[WZ], 2);

	doublereal WW[6];
	WW[VX] = copysign(std::sqrt(Vxp2), W[VX]);
	WW[VY] = -W[VX]*tan_y1 - d34*W[WZ];
	WW[VZ] = W[VZ];
	WW[WX] = W[WX];
	WW[WY] = W[WY];
	WW[WZ] = W[WZ];

	c81_aerod2_u(WW, &VAM, TNG, &OUTA, const_cast<c81_data *>(data), unsteadyflag);
	doublereal cx_1, cy_1, cmz_1;
	if (qD > std::numeric_limits<doublereal>::epsilon()) {
		cx_1 = TNG[0]/(qD*chord);
		cy_1 = TNG[1]/(qD*chord);
		cmz_1 = TNG[5]/(qD*chord*chord);
	}else{
		cx_1 = 0.;
		cy_1 = 0.;
		cmz_1 = 0.;
	}

	doublereal dCd_alpha = (cx_1-cx_0)/dDeltaY1;
	doublereal dCl_alpha = (cy_1-cy_0)/dDeltaY1;
	doublereal dCm_alpha = (cmz_1-cmz_0)/dDeltaY1;

	doublereal C11 = (A1+A2)*b1*b2*d*d;
	doublereal C12 = (A1*b1+A2*b2)*d;

	doublereal qDc = qD*chord;

	fq.Put(1, 1, qDc*dCd_alpha*C11);
	fq.Put(1, 2, qDc*dCd_alpha*C12);
	fq.Put(2, 1, qDc*dCl_alpha*C11);
	fq.Put(2, 2, qDc*dCl_alpha*C12);

	cq.Put(3, 1, qDc*chord*dCm_alpha*C11);
	cq.Put(3, 2, qDc*chord*dCm_alpha*C12);	

	/* computing the J matrix */
	/* f_{/\tilde{v}} */
	doublereal dY_V_11, dY_V_12;

	/* (C_{/Uinf}*q + D_{/Uinf}*u)*Uinf_{v} */
	dY_V_11 = (((A1+A2)*b1*b2*d*(4./chord)*q1 + (A1*b1+A2*b2)*(2./chord)*q2)*W[VX]/Uinf) + ((1-A1-A2)*dU_V_11);
	dY_V_12 = (((A1+A2)*b1*b2*d*(4./chord)*q1 + (A1*b1+A2*b2)*(2./chord)*q2)*W[VY]/Uinf) + ((1-A1-A2)*dU_V_12);

	doublereal dCfy_Uinf = 0.5*clalpha*(-dot_alpha_pivot + (chord*a*ddot_alpha)/Uinf)/(d*Uinf);
	doublereal rho = VAM.density;

	doublereal cy = cy_0 + clalpha/2/d*( dot_alpha_pivot - a/d*ddot_alpha);

	J.Put(1, 1, rho*chord*cx_0*W[VX] + qDc*( dCd_alpha*dY_V_11));
	J.Put(1, 2, rho*chord*cx_0*W[VY] + qDc*( dCd_alpha*dY_V_12));
	J.Put(1, 3, rho*chord*cx_0*W[VZ]);
	J.Put(2, 1, rho*chord*cy*W[VX] + qDc*( dCl_alpha*dY_V_11 + dCfy_Uinf*W[VX]/Uinf) );
	J.Put(2, 2, rho*chord*cy*W[VY] + qDc*( dCl_alpha*dY_V_12 + dCfy_Uinf*W[VY]/Uinf ) );
	J.Put(2, 3, rho*chord*cy*W[VZ]);
	J.Put(3, 1, rho*chord*cz_0*W[VX]);
	J.Put(3, 1, rho*chord*cz_0*W[VY]);
	J.Put(3, 3, rho*chord*cz_0*W[VZ]);

	/* f_{/\tilde{omega}} */
	J.Put(1, 6, qDc*dCd_alpha*(1.-A1-A2)*dU_W_13);
	J.Put(2, 6, qDc*dCl_alpha*(1.-A1-A2)*dU_W_13);

	/* c_{/\tilde{v}} */
	doublereal dCmz_Uinf = 0.5*clalpha*(dot_alpha_pivot - ( ((chord*a)/(Uinf)) - ((chord)/(4*Uinf)) )*ddot_alpha +dot_alpha)/(d*4*Uinf);

	doublereal cmz = cmz_0 + clalpha/2/d*(-dot_alpha_pivot/4 + (a/4/d - 1/d/16)*ddot_alpha - dot_alpha/4);

	J.Put(4, 1, rho*chord*chord*cmx_0*W[VX]);
	J.Put(4, 2, rho*chord*chord*cmx_0*W[VY]);
	J.Put(4, 3, rho*chord*chord*cmx_0*W[VZ]);
	J.Put(5, 1, rho*chord*chord*cmy_0*W[VX]);
	J.Put(5, 2, rho*chord*chord*cmy_0*W[VY]);
	J.Put(5, 3, rho*chord*chord*cmy_0*W[VZ]);
	J.Put(6, 1, rho*chord*chord*cmz*W[VX] +qDc*chord*( dCm_alpha*dY_V_11 + dCmz_Uinf*W[VX]/Uinf ) );
	J.Put(6, 2, rho*chord*chord*cmz*W[VY] +qDc*chord*( dCm_alpha*dY_V_12 + dCmz_Uinf*W[VY]/Uinf ) );
	J.Put(6, 3, rho*chord*chord*cmz*W[VZ]);

	/* c_{/\tilde{omega}} */
	J.Put(6, 6, qDc*chord*dCm_alpha*(1.-A1-A2)*dU_W_13 );

#ifdef  DEBUG_JACOBIAN_UNSTEADY	

	printf("G/v matrix\n");
	printf("%lf %lf %lf\n", 0., 0., 0.);
	printf("%lf %lf %lf\n", dG_V_21, dG_V_22, 0.);

	printf("G/w matrix\n");
	printf("%lf %lf %lf\n", 0., 0., 0.);
	printf("%lf %lf %lf\n", 0., 0., dG_W_23);

	printf("f/q matrix\n");
	printf("%lf %lf \n", fq(1,1), fq(1,2));
	printf("%lf %lf \n", fq(2,1), fq(2,2));
	printf("%lf %lf \n", 0., 0.);

	printf("c/q matrix\n");
	printf("%lf %lf\n", 0., 0.);
	printf("%lf %lf\n", 0., 0.);
	printf("%lf %lf\n", cq(3,1), cq(3,2));

	printf("J matrix\n");
	doublereal iii, jjj;
	for( iii=1; iii<=6; iii++){
		for( jjj=1; jjj<=6; jjj++){
			printf("%lf ", J(iii,jjj));
		}
		printf("\n");
	}
	
	FILE *fd;
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/X.mat","w");
	printf("DATI X\n");
	for( iii=1; iii<=2; iii++){
		printf("%lf ", XCurr(iFirstIndex+iii));
		fprintf(fd,"%15.7e ", XCurr(iFirstIndex+iii));
	}
	fclose(fd);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/XP.mat","w");
	printf("\n");
	printf("DATI XP\n");
	for( iii=1; iii<=2; iii++){
		printf("%lf ", XPrimeCurr(iFirstIndex+iii));
		fprintf(fd,"%15.7e ", XPrimeCurr(iFirstIndex+iii));
	}
	printf("\n");
	fclose(fd);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/W.mat","w");
	printf("DATI W\n");
	printf("%lf %lf %lf %lf %lf %lf", W[VX], W[VY], W[VZ], W[WX], W[WY], W[WZ]);
	fprintf(fd,"%15.7e %15.7e %15.7e %15.7e %15.7e %15.7e", W[VX], W[VY], W[VZ], W[WX], W[WY], W[WZ]);
	fclose(fd);
	printf("\n");
	printf("PARAMETRI \n");
	printf("d14 %lf ", d14);
	printf("\nd34 %lf ", d34);
	printf("\nA1 %lf ", A1);
	printf("\nA2 %lf ", A2);
	printf("\nb1 %lf ", b1);
	printf("\nb2 %lf ", b2);
	printf("\nchord %lf ", chord);
	printf("\na %lf ", a);
	printf("\nrho %lf ", rho);
	printf("\nomegaPD %lf ", omegaPD);
	printf("\ncx_0 %lf ", cx_0);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/cx_0.mat","w");
	fprintf(fd,"%15.7e", cx_0);
	fclose(fd);
	printf("\ncy_0 %lf ", cy_0);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/cy_0.mat","w");
	fprintf(fd,"%15.7e ", cy_0);
	fclose(fd);
	printf("\ncz_0 %lf ", cz_0);
	printf("\ncmx_0 %lf ", cmx_0);
	printf("\ncmy_0 %lf ", cmy_0);
	printf("\ncmz_0 %lf ", cmz_0);
	printf("\ndCl_alpha %lf ", dCl_alpha);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/dCl_alpha.mat","w");
	fprintf(fd,"%15.7e", dCl_alpha);
	fclose(fd);
	printf("\ndCd_alpha %lf", dCd_alpha);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/dCd_alpha.mat","w");
	fprintf(fd,"%15.7e", dCd_alpha);
	fclose(fd);
	printf("\ndCm_alpha %lf", dCm_alpha);
	printf("\nclalpha %lf ", clalpha);
	printf("\n q1 q2 i %lf %lf %d", q1, q2, i);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/ddot_alpha.mat","w");
	fprintf(fd,"%15.7e", ddot_alpha);
	printf("%lf", ddot_alpha);
	fclose(fd);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/dot_alpha.mat","w");
	fprintf(fd,"%15.7e", dot_alpha);
	printf("%lf", dot_alpha);
	fclose(fd);
	fd = fopen("/home/mattia/work/mbdyn/mbdyn/aero/UnsteadyAeroJ_test/dot_alpha_pivot.mat","w");
	fprintf(fd,"%15.7e", dot_alpha_pivot);
	printf("%lf", dot_alpha_pivot);
	fclose(fd);
	getchar();	
#endif /* DEBUG_JACOBIAN_UNSTEADY */	

	AeroData::GetForcesJacForwardDiff_int(i, W, TNG, J, OUTA);

	// probably, we need to reset fq, cq
}
void
C81TheodorsenAeroData::AfterConvergence(int i,	
	const VectorHandler& X, const VectorHandler& XP )
{
	/* aggiorno i valori delle variabili per il calcolo
 	 * della derivata attraverso le differenze finite */
	prev_alpha_pivot = alpha_pivot;
	prev_dot_alpha = dot_alpha;
	prev_time = pTime->dGet();

}

/* C81TheodorsenAeroData - end */

#endif /* USE_UNSTEADY_2STATES */
