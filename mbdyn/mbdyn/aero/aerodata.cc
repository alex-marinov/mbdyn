/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <mynewmem.h>
#include <aerodata.h>
extern "C" {
#include <aerod2.h>
};

AeroMemory::AeroMemory(DriveCaller *pt)
: a(0), t(0), iPoints(0), pTime(pt), numUpdates(0)
{
	NO_OP;
}

AeroMemory::~AeroMemory(void)
{
	if (iPoints > 0) {
		ASSERT(pTime);
		ASSERT(a);

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
			THROW(ErrGeneric());
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
   	VAM[0] = rho;
   	VAM[1] = c;
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
		THROW(ErrGeneric());
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
   	VAM[2] = chord;
   	VAM[3] = forcepoint;
   	VAM[4] = velocitypoint;
   	VAM[5] = twist;
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
      		THROW(ErrGeneric());
   	}
	
	return RestartUnsteady(out);
}
   
int 
STAHRAeroData::GetForces(int i, doublereal* W, doublereal* TNG, 
		doublereal* OUTA) 
{
	switch (unsteadyflag) {
	case AeroData::HARRIS: 
	case AeroData::BIELAWA: {
		Predict(i, atan2(-W[1], W[0]), OUTA[ALF1], OUTA[ALF2]);
		break;
	}

	default:
		break;
	}

	integer u = unsteadyflag;
   	__FC_DECL__(aerod2)(W, VAM, TNG, OUTA, &u, &Omega, &profile);
   	return 0;
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
C81AeroData::GetForces(int i, doublereal* W, doublereal* TNG, doublereal* OUTA) 
{
	switch (unsteadyflag) {
	case AeroData::HARRIS: 
	case AeroData::BIELAWA: {
		Predict(i, atan2(-W[1], W[0]), OUTA[ALF1], OUTA[ALF2]);
		break;
	}

	default:
		break;
	}

   	return c81_aerod2_u(W, VAM, TNG, OUTA, (c81_data*)data, unsteadyflag);
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

	for (int i = nprofiles-1; i--; ) {
		if (abscissa > upper_bounds[i]) {
			curr_data = i+1;
			return;
		}
	}

	curr_data = 0;
}

int
C81MultipleAeroData::GetForces(int i, doublereal* W, doublereal* TNG, doublereal* OUTA) 
{
	switch (unsteadyflag) {
	case AeroData::HARRIS: 
	case AeroData::BIELAWA: {
		Predict(i, atan2(-W[1], W[0]), OUTA[ALF1], OUTA[ALF2]);
		break;
	}

	default:
		break;
	}

   	return c81_aerod2_u(W, VAM, TNG, OUTA, (c81_data*)data[curr_data],
			unsteadyflag);
}


