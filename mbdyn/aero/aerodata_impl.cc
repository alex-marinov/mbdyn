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
#ifdef USE_AEROD2_F
#include "aerod2.h"
#endif // USE_AEROD2_F
#include "c81data.h"

#ifdef USE_AEROD2_F

/* STAHRAeroData - begin */

STAHRAeroData::STAHRAeroData(int i_p, int i_dim,
	AeroData::UnsteadyModel u, integer p,
	DriveCaller *ptime)
: AeroData(i_p, i_dim, u, ptime),
profile(p)
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

#endif // USE_AEROD2_F

/* C81AeroData - begin */

C81AeroData::C81AeroData(int i_p, int i_dim,
	AeroData::UnsteadyModel u, integer p,
	const c81_data* d, DriveCaller *ptime)
: AeroData(i_p, i_dim, u, ptime),
profile(p), data(d)
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
	int i_p, int i_dim,
	AeroData::UnsteadyModel u,
	std::vector<unsigned>& p,
	std::vector<doublereal>& ub,
	std::vector<const c81_data *>& d,
	DriveCaller *ptime)
: AeroData(i_p, i_dim, u, ptime),
profiles(p), upper_bounds(ub), data(d)
{
	ASSERT(!profiles.empty());
	ASSERT(!upper_bounds.empty());
	ASSERT(!data.empty());
	ASSERT(profiles.size() == upper_bounds.size());
	ASSERT(profiles.size() == data.size());
}

C81MultipleAeroData::~C81MultipleAeroData(void)
{
	NO_OP;
}

std::ostream&
C81MultipleAeroData::Restart(std::ostream& out) const
{
	out << "C81, multiple";
	for (unsigned i = 0; i < profiles.size(); i++) {
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
	const doublereal& omega)
{
	ASSERT(abscissa >= -1. && abscissa <= 1.);

	AeroData::SetSectionData(abscissa, chord, forcepoint, velocitypoint,
		twist, omega);

	for (int i = profiles.size() - 1; i--; ) {
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

	ASSERT(i >= 0);
	ASSERT(unsigned(i) < data.size());

	return c81_aerod2_u(const_cast<doublereal *>(W), &VAM, TNG, &OUTA,
		const_cast<c81_data *>(data[curr_data]), unsteadyflag);
}

int
C81MultipleAeroData::GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
	return AeroData::GetForcesJacForwardDiff_int(i, W, TNG, J, OUTA);
}

/* C81MultipleAeroData - end */

/* C81InterpolatedAeroData - begin */

C81InterpolatedAeroData::C81InterpolatedAeroData(
		integer i_p, int i_dim,
		AeroData::UnsteadyModel u,
		std::vector<unsigned>& p,
		std::vector<doublereal>& ub,
		std::vector<const c81_data *>& d,
		doublereal dcltol,
		DriveCaller *ptime
)
: AeroData(i_p, i_dim, u, ptime),
profiles(p), upper_bounds(ub), data(d),
i_data(i_p*i_dim)
{
	ASSERT(!profiles.empty());
	ASSERT(!upper_bounds.empty());
	ASSERT(!data.empty());
	ASSERT(profiles.size() == upper_bounds.size());
	ASSERT(profiles.size() == data.size());

	ASSERT(i_dim >= 1);
	ASSERT(i_dim <= 3);
	ASSERT(!i_data.empty());

	const doublereal bs[3][4] = {
		{ -1., 1. },
		{ -1., 0., 1. },
		{ -1., -1./std::sqrt(3.), 1./std::sqrt(3.), 1. }
	};

	int i_point = 0;
	for (int b = 0; b < i_dim; b++) {
		GaussDataIterator GDI(i_p);
		PntWght PW = GDI.GetFirst();
		do {
			doublereal dCsi = PW.dGetPnt();

			doublereal bl = bs[i_dim - 1][b + 1] - bs[i_dim - 1][b];
			doublereal bm = bs[i_dim - 1][b + 1] + bs[i_dim - 1][b];
			dCsi = (bm + dCsi*bl)/2.;

			if (c81_data_merge(data.size(), &data[0], &upper_bounds[0],
				dCsi, dcltol, &i_data[i_point]))
			{
				// function logs error
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			i_point++;
		} while (GDI.fGetNext(PW));
	}
}

C81InterpolatedAeroData::~C81InterpolatedAeroData(void)
{
	for (unsigned i = 0; i < i_data.size(); i++) {
		c81_data_destroy(&i_data[i]);
	}
}

std::ostream&
C81InterpolatedAeroData::Restart(std::ostream& out) const
{
	out << "C81, interpolated";
	for (unsigned i = 0; i < profiles.size(); i++) {
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
	const doublereal& omega)
{
	ASSERT(abscissa >= -1. && abscissa <= 1.);

	AeroData::SetSectionData(abscissa, chord, forcepoint, velocitypoint,
		twist, omega);
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
		const_cast<c81_data *>(&i_data[i]), unsteadyflag);
}

int
C81InterpolatedAeroData::GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
	return AeroData::GetForcesJacForwardDiff_int(i, W, TNG, J, OUTA);
}

/* C81InterpolatedAeroData - end */

/* TheodorsenAeroData - begin */

static const doublereal TheodorsenParams[2][4] = {
	{ 0.165, 0.335, 0.0455, 0.3 },
	{ 0.165, 0.335, 0.041, 0.32 }
};

TheodorsenAeroData::TheodorsenAeroData(
	int i_p, int i_dim,
	AeroData *pa,
	DriveCaller *ptime)
: AeroData(i_p, i_dim, STEADY, ptime),
iParam(0),
d14(0.), d34(0.),
chord(-1.),
a(0.),
A1(0.), A2(0.), b1(0.), b2(0.),
alpha_pivot(0), dot_alpha_pivot(0), dot_alpha(0), ddot_alpha(0),
cfx_0(0), cfy_0(0), cmz_0(0), clalpha(0),
prev_alpha_pivot(0), prev_dot_alpha(0),
prev_time(pTime->dGet()),
pAeroData(pa)
{
	ASSERT(pAeroData->Unsteady() == STEADY);
	ASSERT(pAeroData->iGetNumDof() == 0);
}

TheodorsenAeroData::~TheodorsenAeroData(void)
{
	if (alpha_pivot) {
		SAFEDELETEARR(alpha_pivot);
	}

	if (pAeroData) {
		SAFEDELETE(pAeroData);
	}
}

std::ostream&
TheodorsenAeroData::Restart(std::ostream& out) const
{
	out << "theodorsen";

	return pAeroData->Restart(out);
}

void
TheodorsenAeroData::SetAirData(const doublereal& rho, const doublereal& c)
{
	AeroData::SetAirData(rho, c);
	pAeroData->SetAirData(rho, c);
}

void
TheodorsenAeroData::SetSectionData(const doublereal& abscissa,
	const doublereal& chord,
	const doublereal& forcepoint,
	const doublereal& velocitypoint,
	const doublereal& twist,
	const doublereal& omega)
{
	AeroData::SetSectionData(abscissa, chord, forcepoint, velocitypoint, twist, omega);
	pAeroData->SetSectionData(abscissa, chord, forcepoint, velocitypoint, twist, omega);
}

// aerodynamic models with internal states
unsigned int
TheodorsenAeroData::iGetNumDof(void) const
{
	return 2;
}

DofOrder::Order
TheodorsenAeroData::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

void
TheodorsenAeroData::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr,
	integer iFirstIndex, integer iFirstSubIndex,
	int i, const doublereal* W, doublereal* TNG, outa_t& OUTA)
{
	// FIXME: should do it earlier
	if (alpha_pivot == 0) {
		doublereal *d = 0;
		SAFENEWARR(d, doublereal, 10*GetNumPoints());

#ifdef HAVE_MEMSET
		memset(d, '0', 10*GetNumPoints()*sizeof(doublereal));
#else // ! HAVE_MEMSET
		for (int i = 0; i < 10*GetNumPoints(); i++) {
			d[i] = 0.;
		}
#endif // ! HAVE_MEMSET

		alpha_pivot = d;
		d += GetNumPoints();
		dot_alpha_pivot = d;
		d += GetNumPoints();
		dot_alpha = d;
		d += GetNumPoints();
		ddot_alpha = d;
		d += GetNumPoints();
		cfx_0 = d;
		d += GetNumPoints();
		cfy_0 = d;
		d += GetNumPoints();
		cmz_0 = d;
		d += GetNumPoints();
		clalpha = d;
		d += GetNumPoints();
		prev_alpha_pivot = d;
		d += GetNumPoints();
		prev_dot_alpha = d;
		d += GetNumPoints();
	}

	doublereal q1 = XCurr(iFirstIndex + 1);
	doublereal q2 = XCurr(iFirstIndex + 2);

	doublereal q1p = XPrimeCurr(iFirstIndex + 1);
	doublereal q2p = XPrimeCurr(iFirstIndex + 2);

	d14 = VAM.force_position;
	d34 = VAM.bc_position;
	if (std::abs(d34 - d14) < std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("TheodorsenAeroData::AssRes() [#" << i << "]: aerodynamic center and boundary condition point are almost coincident" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	chord = VAM.chord;

	a = (d14 + d34)/chord;

	//doublereal Uinf = sqrt(W[VX]*W[VX] + W[VY]*W[VY]);
	doublereal Uinf = sqrt(W[VX]*W[VX] + (W[VY]+d34*W[WZ])*(W[VY]+d34*W[WZ]));

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
	//WW[VY] = -W[VX]*tan_y1 - d34*W[WZ];
	WW[VY] = -WW[VX]*tan_y1 - d34*W[WZ];
	WW[VZ] = W[VZ];
	WW[WX] = W[WX];
	WW[WY] = W[WY];
	WW[WZ] = W[WZ];

	pAeroData->GetForces(i, WW, TNG, OUTA);

	doublereal UUinf2 = Uinf*Uinf + W[VZ]*W[VZ];
	doublereal qD = .5*VAM.density*UUinf2;
	if (qD > std::numeric_limits<doublereal>::epsilon()) {
		cfx_0[i] = OUTA.cd;
		cfy_0[i] = OUTA.cl;
#if 0
		cfz_0 = 0.;
		cmx_0 = 0.;
		cmy_0 = 0.;
#endif
		cmz_0[i] = OUTA.cm;
	} else {
		cfx_0[i] = 0.;
		cfy_0[i] = 0.;
#if 0
		cfz_0 = 0.;
		cmx_0 = 0.;
		cmy_0 = 0.;
#endif
		cmz_0[i] = 0.;
	}

	alpha_pivot[i] = (u1*d34 - u2*d14)/(d34 - d14);
	dot_alpha[i] = Uinf*((u1 - u2)/(d34 - d14));

	doublereal Delta_t = pTime->dGet() - prev_time;
	/* controllo che dovrebbe servire solo al primo istante di
 	 * di tempo, in cui non calcolo le derivate per differenze
 	 * finite per non avere 0 a denominatore. */ 
	if (Delta_t > std::numeric_limits<doublereal>::epsilon()) {
		dot_alpha_pivot[i] = (alpha_pivot[i] - prev_alpha_pivot[i])/Delta_t;
		ddot_alpha[i] = (dot_alpha[i] - prev_dot_alpha[i])/Delta_t;

	} else {
		dot_alpha_pivot[0] = 0.;
		ddot_alpha[0] = 0.;
	}

	clalpha[i] = OUTA.clalpha;
	if (clalpha[i] > 0.) {
		//if (Uinf > std::numeric_limits<doublereal>::epsilon()) {		
		if (Uinf > 1.e-5) {		
			doublereal cl = OUTA.cl + clalpha[i]/2/d*(dot_alpha_pivot[i] - a/d*ddot_alpha[i]);
			doublereal cd = OUTA.cd;
			doublereal cm = OUTA.cm + clalpha[i]/2/d*(-dot_alpha_pivot[i]/4 + (a/4/d - 1/d/16)*ddot_alpha[i] - dot_alpha[i]/4);
			doublereal v[3], vp;
			v[0] = W[0];
			v[1] = W[1] + d34*W[5];
			v[2] = W[2] - d34*W[4];
			vp = sqrt(v[0]*v[0] + v[1]*v[1]);
			TNG[0] = -qD*chord*(cl*v[1] + cd*v[0])/vp;
       			TNG[1] = qD*chord*(cl*v[0] - cd*v[1])/vp;
			TNG[5] = qD*chord*chord*cm + d14*TNG[1];
		} else {
			TNG[0] = 0.;
			TNG[1] = 0.;
			TNG[2] = 0.;
			TNG[3] = 0.;
			TNG[4] = 0.;
			TNG[5] = 0.;
		}
	} else {
		TNG[0] = 0.;
		TNG[1] = 0.;
		TNG[2] = 0.;
		TNG[3] = 0.;
		TNG[4] = 0.;
		TNG[5] = 0.;
	}

	WorkVec.PutCoef(iFirstSubIndex + 1, -q1p + q2);
	WorkVec.PutCoef(iFirstSubIndex + 2, -q2p - b1*b2*d*d*q1 - (b1 + b2)*d*q2 + u2);
}

void
TheodorsenAeroData::AssJac(FullSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr,
	integer iFirstIndex, integer iFirstSubIndex,
	const Mat3xN& vx, const Mat3xN& wx, Mat3xN& fq, Mat3xN& cq,
	int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA)
{
#ifdef DEBUG_JACOBIAN_UNSTEADY
	doublereal q1 = XCurr(iFirstIndex + 1);
	doublereal q2 = XCurr(iFirstIndex + 2);
#endif // DEBUG_JACOBIAN_UNSTEADY

	doublereal Uinf = sqrt(W[VX]*W[VX] + W[VY]*W[VY]);
	doublereal d = 2*Uinf/chord;
#ifdef DEBUG_JACOBIAN_UNSTEADY
	doublereal UUinf2 = Uinf*Uinf + W[VZ]*W[VZ];
#endif // DEBUG_JACOBIAN_UNSTEADY
#if 0
	doublereal qD = .5*VAM.density*UUinf2;
#endif

	// doublereal u1 = atan2(- W[VY] - W[WZ]*d14, W[VX]);
#ifdef DEBUG_JACOBIAN_UNSTEADY
	doublereal u2 = atan2(- W[VY] - W[WZ]*d34, W[VX]);
#endif // DEBUG_JACOBIAN_UNSTEADY

#if 0
	doublereal y1 = (A1 + A2)*b1*b2*d*d*q1 + (A1*b1 + A2*b2)*d*q2
		+ (1 - A1 - A2)*u2;
#endif

	WorkMat.IncCoef(iFirstSubIndex + 1, iFirstSubIndex + 1, 1.);
	WorkMat.IncCoef(iFirstSubIndex + 2, iFirstSubIndex + 2, 1.);

	WorkMat.IncCoef(iFirstSubIndex + 1, iFirstSubIndex + 2, -dCoef);
	WorkMat.IncCoef(iFirstSubIndex + 2, iFirstSubIndex + 1, dCoef*b1*b2*d*d);
	WorkMat.IncCoef(iFirstSubIndex + 2, iFirstSubIndex + 2, dCoef*(b1 + b2)*d);

#if 0
	//if (Uinf > std::numeric_limits<doublereal>::epsilon()) {		
	if (Uinf > 1.e-2) {		
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
		int iIndexColumn;
		/* faccio i calcoli tenendo conto della sparsità della matrice g_{\/\tilde{v}} */
		#if 1
		for (iIndexColumn = 1; iIndexColumn <= vx.iGetNumCols(); iIndexColumn++){
			WorkMat.IncCoef(iFirstSubIndex+2, iIndexColumn, dG_V_21*vx(1,iIndexColumn) + dG_V_22*vx(2,iIndexColumn));
			WorkMat.IncCoef(iFirstSubIndex+2, iIndexColumn, dG_W_23*wx(3,iIndexColumn));
		}
		#endif
		/* computing the matrix fq (f_a_{/q}) 3x2 */
	
		/* computing the derivative of the aerodynamic coefficient in the lookup table
 		 * using the finite difference method */
		/* perturbazione per calcolo delle derivare con le differenze finite */
		doublereal dDeltaY1 = 1.*M_PI/180.;
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
			cx_1 = OUTA.cd;
			cy_1 = OUTA.cl;
			cmz_1 = OUTA.cm;
		}else{
			cx_1 = 0.;
			cy_1 = 0.;
			cmz_1 = 0.;
		}
	
		doublereal dCd_alpha = (cx_1-cfx_0[i])/dDeltaY1;
		doublereal dCl_alpha = (cy_1-cfy_0[i])/dDeltaY1;
		doublereal dCm_alpha = (cmz_1-cmz_0[i])/dDeltaY1;
	
		doublereal C11 = (A1+A2)*b1*b2*d*d;
		doublereal C12 = (A1*b1+A2*b2)*d;
	
		doublereal qDc = qD*chord;
		/* ATTENZIONE: le derivate aerodinamiche calcolate (vedi tecman)
 		 * sono le dirivate di L D e M34 e vanno quindi routate per avere
 		 * le derivate di TNG */
		doublereal sin_alpha = (W[VY]+d34*W[WZ])/sqrt(W[VX]*W[VX]+W[VY]*W[VY]);
		doublereal cos_alpha = (W[VX])/sqrt(W[VX]*W[VX]+W[VY]*W[VY]);
	
		doublereal D_q1 = qDc*dCd_alpha*C11;
		doublereal D_q2 = qDc*dCd_alpha*C12;
		doublereal L_q1 = qDc*dCl_alpha*C11;
		doublereal L_q2 = qDc*dCl_alpha*C12;
		
		doublereal M_q1 = qDc*chord*dCm_alpha*C11;
		doublereal M_q2 = qDc*chord*dCm_alpha*C12;
		#if 0
		fq.Put(1, 1, (-L_q1*sin_alpha -D_q1*cos_alpha));
		fq.Put(1, 2, (-L_q2*sin_alpha -D_q2*cos_alpha));
		fq.Put(2, 1, (L_q1*cos_alpha -D_q1*sin_alpha));
		fq.Put(2, 2, (L_q2*cos_alpha -D_q2*sin_alpha));
		
		cq.Put(3, 1, (M_q1 +d14*fq(2,1)));
		cq.Put(3, 2, (M_q2 +d14*fq(2,2)));
		#endif
		/* computing the J matrix */
		/* f_{/\tilde{v}} */
		doublereal dY_V_11, dY_V_12;
	
		/* (C_{/Uinf}*q + D_{/Uinf}*u)*Uinf_{v} */
		dY_V_11 = (((A1+A2)*b1*b2*d*(4./chord)*q1 + (A1*b1+A2*b2)*(2./chord)*q2)*W[VX]/Uinf) + ((1-A1-A2)*dU_V_11);
		dY_V_12 = (((A1+A2)*b1*b2*d*(4./chord)*q1 + (A1*b1+A2*b2)*(2./chord)*q2)*W[VY]/Uinf) + ((1-A1-A2)*dU_V_12);
	
		doublereal dCfy_Uinf = 0.5*clalpha[i]*(-dot_alpha_pivot[i] + (chord*a*ddot_alpha[i])/Uinf)/(d*Uinf);
		doublereal rho = VAM.density;
	
		doublereal cy = cfy_0[i] + clalpha[i]/2/d*( dot_alpha_pivot[i] - a/d*ddot_alpha[i]);
	
		Mat6x6 Jaero = 0.;
		Jaero(1,1) = rho*chord*cfx_0[i]*W[VX] + qDc*( dCd_alpha*dY_V_11);
		Jaero(1,2) = rho*chord*cfx_0[i]*W[VY] + qDc*( dCd_alpha*dY_V_12);
		Jaero(1,3) = rho*chord*cfx_0[i]*W[VZ];
		Jaero(2,1) = rho*chord*cy*W[VX] + qDc*( dCl_alpha*dY_V_11 + dCfy_Uinf*W[VX]/Uinf);
		Jaero(2,2) = rho*chord*cy*W[VY] + qDc*( dCl_alpha*dY_V_12 + dCfy_Uinf*W[VY]/Uinf );
		Jaero(2,3) = rho*chord*cy*W[VZ];
		Jaero(3,1) = rho*chord*cfz_0*W[VX];
		Jaero(3,1) = rho*chord*cfz_0*W[VY];
		Jaero(3,3) = rho*chord*cfz_0*W[VZ];
	
		/* f_{/\tilde{omega}} */
		Jaero(1,6) = qDc*dCd_alpha*(1.-A1-A2)*dU_W_13;
		Jaero(2,6) = qDc*dCl_alpha*(1.-A1-A2)*dU_W_13;
	
		/* c_{/\tilde{v}} */
		doublereal dCmz_Uinf = 0.5*clalpha[i]*(dot_alpha_pivot[i] - ( ((chord*a)/(Uinf)) - ((chord)/(4*Uinf)) )*ddot_alpha[i] + dot_alpha[i])/(d*4*Uinf);
	
		doublereal cmz = cmz_0 + clalpha[i]/2/d*(-dot_alpha_pivot[i]/4 + (a/4/d - 1/d/16)*ddot_alpha[i] - dot_alpha[i]/4);
	
		Jaero(4,1) = rho*chord*chord*cmx_0*W[VX];
		Jaero(4,2) = rho*chord*chord*cmx_0*W[VY];
		Jaero(4,3) = rho*chord*chord*cmx_0*W[VZ];
		Jaero(5,1) = rho*chord*chord*cmy_0*W[VX];
		Jaero(5,2) = rho*chord*chord*cmy_0*W[VY];
		Jaero(5,3) = rho*chord*chord*cmy_0*W[VZ];
		Jaero(6,1) = rho*chord*chord*cmz*W[VX] +qDc*chord*( dCm_alpha*dY_V_11 + dCmz_Uinf*W[VX]/Uinf );
		Jaero(6,2) = rho*chord*chord*cmz*W[VY] +qDc*chord*( dCm_alpha*dY_V_12 + dCmz_Uinf*W[VY]/Uinf );
		Jaero(6,3) = rho*chord*chord*cmz*W[VZ];
	
		/* c_{/\tilde{omega}} */
		Jaero(6,6) = qDc*chord*dCm_alpha*(1.-A1-A2)*dU_W_13;
	
		/* Jaero è lo jacobiano delle forze aerodinamiche L, D e M34
		 * occorre ruotarlo per avere lo jacobiano corretto J.
		 * Solo le righe 1 2 e 6 sono coivolte nella rotazione */
		J = Jaero;
		/* calcolo le forze aerodinamiche */
		doublereal Drag = qDc*cfx_0[i];
		doublereal Lift = qDc*cfy_0[i] + qDc*clalpha[i]/2/d*(dot_alpha_pivot[i] - a/d*ddot_alpha[i]);
		/* calcolo le derivate del cos e sen dell'angolo di rotazione */
		doublereal den = (W[VX]*W[VX]+W[VY]*W[VY])*sqrt((W[VX]*W[VX]+W[VY]*W[VY]));
		doublereal sin_alpha_vx = -W[VX]*(W[VY]+d34*W[WZ])/den;
		doublereal sin_alpha_vy = ( W[VX]*W[VX] - d34*W[WZ]*W[VY] )/den;
		doublereal sin_alpha_wz = d34/sqrt(W[VX]*W[VX]+W[VY]*W[VY]);
		doublereal cos_alpha_vx = ( W[VY]*W[VY]  )/den;
		doublereal cos_alpha_vy = -( W[VX]*W[VY]  )/den;
		#if 1
		J(1,1) = -Jaero(2,1)*sin_alpha -Jaero(1,1)*cos_alpha - Lift*sin_alpha_vx - Drag*cos_alpha_vx;
		J(1,2) = -Jaero(2,2)*sin_alpha -Jaero(1,2)*cos_alpha - Lift*sin_alpha_vy - Drag*cos_alpha_vy;
		J(1,3) = -Jaero(2,3)*sin_alpha -Jaero(1,3)*cos_alpha;
		J(1,6) = -Jaero(2,6)*sin_alpha -Jaero(1,6)*cos_alpha - Lift*sin_alpha_wz;
	
		J(2,1) = Jaero(2,1)*cos_alpha -Jaero(1,1)*sin_alpha + Lift*cos_alpha_vx - Drag*sin_alpha_vx;
		J(2,2) = Jaero(2,2)*cos_alpha -Jaero(1,2)*sin_alpha + Lift*cos_alpha_vy - Drag*sin_alpha_vy;
		J(2,3) = Jaero(2,3)*cos_alpha -Jaero(1,3)*sin_alpha;
		J(2,6) = Jaero(2,6)*cos_alpha -Jaero(1,6)*sin_alpha - Drag*sin_alpha_wz;
	
		J(6,1) = Jaero(6,1) + d14*J(2,1);
		J(6,2) = Jaero(6,2) + d14*J(2,2);
		J(6,3) = Jaero(6,3) + d14*J(2,3);
		J(6,4) = d14*J(2,4);
		J(6,5) = d14*J(2,5);
		J(6,6) = Jaero(6,6) + d14*J(2,6);
		#endif
	}
#endif
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
	for( int iii=1; iii<=6; iii++){
		for( int jjj=1; jjj<=6; jjj++){
			printf("%lf ", J(iii,jjj));
		}
		printf("\n");
	}
	
	FILE *fd;
	fd = fopen("X.mat","w");
	printf("DATI X\n");
	for( int iii=1; iii<=2; iii++){
		printf("%lf ", XCurr(iFirstIndex+iii));
		fprintf(fd,"%15.7e ", XCurr(iFirstIndex+iii));
	}
	fclose(fd);
	fd = fopen("XP.mat","w");
	printf("\n");
	printf("DATI XP\n");
	for( int iii=1; iii<=2; iii++){
		printf("%lf ", XPrimeCurr(iFirstIndex+iii));
		fprintf(fd,"%15.7e ", XPrimeCurr(iFirstIndex+iii));
	}
	printf("\n");
	fclose(fd);
	fd = fopen("W.mat","w");
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
	printf("\ncfx_0 %lf ", cfx_0[i]);
	fd = fopen("cfx_0.mat","w");
	fprintf(fd,"%15.7e", cfx_0[i]);
	fclose(fd);
	printf("\ncfy_0 %lf ", cfy_0[i]);
	fd = fopen("cfy_0.mat","w");
	fprintf(fd,"%15.7e ", cfy_0[i]);
	fclose(fd);
	printf("\ncmz_0 %lf ", cmz_0[i]);
	printf("\ndCl_alpha %lf ", dCl_alpha);
	fd = fopen("dCl_alpha.mat","w");
	fprintf(fd,"%15.7e", dCl_alpha);
	fclose(fd);
	printf("\ndCd_alpha %lf", dCd_alpha);
	fd = fopen("dCd_alpha.mat","w");
	fprintf(fd,"%15.7e", dCd_alpha);
	fclose(fd);
	printf("\ndCm_alpha %lf", dCm_alpha);
	printf("\nclalpha %lf ", clalpha[i]);
	printf("\n q1 q2 i %lf %lf %d", q1, q2, i);
	fd = fopen("ddot_alpha.mat","w");
	fprintf(fd,"%15.7e", ddot_alpha[i]);
	printf("%lf", ddot_alpha[i]);
	fclose(fd);
	fd = fopen("dot_alpha.mat","w");
	fprintf(fd,"%15.7e", dot_alpha[i]);
	printf("%lf", dot_alpha[i]);
	fclose(fd);
	fd = fopen("dot_alpha_pivot.mat","w");
	fprintf(fd,"%15.7e", dot_alpha_pivot[i]);
	printf("%lf", dot_alpha_pivot[i]);
	fclose(fd);
	getchar();	
#endif /* DEBUG_JACOBIAN_UNSTEADY */	

	pAeroData->GetForcesJac(i, W, TNG, J, OUTA);

	// probably, we need to reset fq, cq
}

void
TheodorsenAeroData::AfterConvergence(int i,	
	const VectorHandler& X, const VectorHandler& XP )
{
	/* aggiorno i valori delle variabili per il calcolo
 	 * della derivata attraverso le differenze finite */
	prev_alpha_pivot[i] = alpha_pivot[i];
	prev_dot_alpha[i] = dot_alpha[i];

	// same for all
	prev_time = pTime->dGet();
}

/* TheodorsenAeroData - end */
