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
#include "aerodata.h"
#include "gauss.h"
#include "submat.h"

#include "aerod2.h"

AeroMemory::AeroMemory(DriveCaller *pt)
: a(0), t(0), iPoints(0), pTime(pt), numUpdates(0)
{
	NO_OP;
}

AeroMemory::~AeroMemory(void)
{
	if (iPoints > 0) {
		ASSERT(pTime != 0);
		ASSERT(a != 0);

		SAFEDELETE(pTime);
		SAFEDELETEARR(a);
	}
}

#define USE_POLCOE
void
AeroMemory::Predict(int i, doublereal alpha, doublereal &alf1, doublereal &alf2)
{
	/* FIXME: this should be s, but I don't want a malloc here */
	doublereal 	coe[3];
	int		s = StorageSize();
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

	if (s > 0) {
		iPoints = i;

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


C81Data::C81Data(unsigned int uLabel)
: WithLabel(uLabel)
{
   	NO_OP;
}

AeroData::AeroData(AeroData::UnsteadyModel u, DriveCaller *ptime)
: AeroMemory(ptime), unsteadyflag(u), Omega(0.)
{
   	if (u != AeroData::STEADY) {
		if (ptime == 0) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
}

AeroData::~AeroData(void)
{
   	NO_OP;
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

C81AeroData::C81AeroData(AeroData::UnsteadyModel u, integer p,
		const c81_data* d, DriveCaller *ptime)
: AeroData(u, ptime), profile(p), data(d)
{
   	ASSERT(data != NULL);
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

#ifdef USE_C81INTERPOLATEDAERODATA

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

	for (int i = nprofiles-1; i--; ) {
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

#endif /* USE_C81INTERPOLATEDAERODATA */

UMDAeroData::UMDAeroData(DriveCaller *ptime)
: AeroData(STEADY, ptime)
{
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

