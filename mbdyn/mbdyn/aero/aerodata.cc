/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

C81Data::C81Data(unsigned int uLabel)
: WithLabel(uLabel) 
{
   	NO_OP;
}

AeroData::AeroData(integer u) 
: unsteadyflag(u), Omega(0.)
{
   	NO_OP;
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
   
STAHRAeroData::STAHRAeroData(integer u, integer p) 
: AeroData(u), profile(p)
{
   	NO_OP;
}
   
ostream& 
STAHRAeroData::Restart(ostream& out) const 
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
	
   	return out << ", unsteady, " << unsteadyflag;
}
   
int 
STAHRAeroData::GetForces(doublereal* W, doublereal* TNG, doublereal* OUTA) 
{
   	__FC_DECL__(aerod2)(W, VAM, TNG, OUTA, 
			    &unsteadyflag, &Omega, &profile);
   	return 0;
}

C81AeroData::C81AeroData(integer u, integer p, const c81_data* d)
: AeroData(u), profile(p), data(d) 
{
   	ASSERT(data != NULL);
}

ostream& 
C81AeroData::Restart(ostream& out) const 
{
   	return out << "C81, " << profile << ", unsteady, " << unsteadyflag;
}

int
C81AeroData::GetForces(doublereal* W, doublereal* TNG, doublereal* OUTA) 
{
   	return c81_aerod2(W, VAM, TNG, OUTA, (c81_data*)data);
}



C81MultipleAeroData::C81MultipleAeroData(
		integer u,
		integer np,
		integer *p,
		doublereal *ub,
		const c81_data** d
)
: AeroData(u), nprofiles(np), profiles(p), upper_bounds(ub), data(d) 
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

ostream& 
C81MultipleAeroData::Restart(ostream& out) const 
{
   	out << "C81, ";
	for (int i = 0; i < nprofiles; i++) {
		out << profiles[i] << ", " << upper_bounds[i] << ", ";
	}
	return out << "unsteady, " << unsteadyflag;
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
	ASSERT(abscissa > 1. || abscissa < -1.);

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
C81MultipleAeroData::GetForces(doublereal* W, doublereal* TNG, doublereal* OUTA) 
{
   	return c81_aerod2(W, VAM, TNG, OUTA, (c81_data*)data[curr_data]);
}


