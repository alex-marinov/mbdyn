/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/* 
 * Copyright 1999-2000 Lamberto Puggelli <puggelli@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#ifndef HFLUID__H
#define HFLUID__H

#include <hfluid.h>

/* IncompressibleHydraulicFluid - begin */

class IncompressibleHydraulicFluid : public HydraulicFluid {
 protected:
   const doublereal dDensity;
   const doublereal dViscosity;
   
 public:
   IncompressibleHydraulicFluid(unsigned int Label, 
				const doublereal& dDensity,
				const doublereal& dViscosity = 0.,
				const doublereal& dPres0 = -1.,
				const doublereal& dTemp0 = -1.)
     : HydraulicFluid(Label, dPres0, dTemp0), 
     dDensity(dDensity), dViscosity(dViscosity) {
	NO_OP;
     };
   IncompressibleHydraulicFluid(const IncompressibleHydraulicFluid& HF)
     : HydraulicFluid(HF), 
     dDensity(HF.dDensity), dViscosity(HF.dViscosity) {
	NO_OP;
     };
   
   ~IncompressibleHydraulicFluid() {
      NO_OP;
   };
   
   /* crea una copia dell'HF */
   HydraulicFluid* pCopy(void) const {
      HydraulicFluid* pHF = NULL;
      SAFENEWWITHCONSTRUCTOR(pHF, 
			     IncompressibleHydraulicFluid,
			     IncompressibleHydraulicFluid(*this));
      return pHF;
   };
   
   doublereal dGetDensity(void) const {
      return dDensity;
   };
   doublereal dGetDensity(const doublereal& /* dPres */ ) const {
      return dDensity;
   };
   doublereal dGetDensity(const doublereal& /* dPres */ , const doublereal& /* dTemp */ ) const {
      return dDensity;
   };
   
   doublereal dGetDensityDPres(void) const {
      return 0.;
   };
   doublereal dGetDensityDPres(const doublereal& /* dPres */ ) const {
      return 0.;
   };
   doublereal dGetDensityDPres(const doublereal& /* dPres */ , const doublereal& /* dTemp */ ) const {
      return 0.;
   };
   
   doublereal dGetDensityDTemp(void) const {
      return 0.;
   };
   doublereal dGetDensityDTemp(const doublereal& /* dPres */ ) const {
      return 0.;
   };
   doublereal dGetDensityDTemp(const doublereal& /* dPres */ , const doublereal& /* dTemp */ ) const {
      return 0.;
   };

   doublereal dGetViscosity(void) const {
      return dViscosity;
   };
   
   doublereal dGetViscosity(const doublereal& /* dPres */ ) const {
      return dViscosity;
   }; 
   
   doublereal dGetViscosity(const doublereal& /* dPres */ , const doublereal& /* dTemp */ ) const {
      return dViscosity;
   };   
};

/* IncompressibleHydraulicFluid - end */


/* LinearCompressibleHydraulicFluid - begin */

class LinearCompressibleHydraulicFluid : public HydraulicFluid {
 protected:
   const doublereal dDensity;
   const doublereal dBeta;
   const doublereal dDensityDPres;
   const doublereal dViscosity;
   
 public:
   LinearCompressibleHydraulicFluid(unsigned int Label, 
				    const doublereal& dDens,
				    const doublereal& dB,
				    const doublereal& dPr0,
				    const doublereal& dVisc = 0.,
				    const doublereal dTp0 = -1.)
     : HydraulicFluid(Label, dPr0, dTp0), 
   dDensity(dDens), 
   dBeta(dB),
   dDensityDPres(0.),
   dViscosity(dVisc) {
      if (dB == 0.) {
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
      (doublereal&)dDensityDPres = dDens/dB;
   };
   
   LinearCompressibleHydraulicFluid(const LinearCompressibleHydraulicFluid& HF)
     : HydraulicFluid(HF), 
   dDensity(HF.dDensity), 
   dBeta(HF.dBeta),
   dDensityDPres(0.),
   dViscosity(HF.dViscosity) {
      if (dBeta == 0.) {
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
      (doublereal&)dDensityDPres = dDensity/dBeta;
   };
   
   virtual ~LinearCompressibleHydraulicFluid() {
      NO_OP;
   };
   
   /* crea una copia dell'HF */
   HydraulicFluid* pCopy(void) const {
      HydraulicFluid* pHF = NULL;
      SAFENEWWITHCONSTRUCTOR(pHF, 
			     LinearCompressibleHydraulicFluid,
			     LinearCompressibleHydraulicFluid(*this));
      return pHF;
   };

   /* densita' */
   virtual doublereal dGetDensity(void) const {
      return dDensity;
   };
   virtual doublereal dGetDensity(const doublereal& dPres) const {
      return dDensity+(dPres-dPres0)*dDensityDPres;
   };
   virtual doublereal dGetDensity(const doublereal& dPres,
			  const doublereal& /* dTemp */) const {
      return dGetDensity(dPres);
   };
   

   virtual doublereal dGetDensityDPres(void) const {
      return dDensityDPres;
   };
   virtual doublereal dGetDensityDPres(const doublereal& /* dPres */ ) const {
      return dDensityDPres;
   };
   virtual doublereal dGetDensityDPres(const doublereal& /* dPres */ , const doublereal& /* dTemp */ ) const {
      return dDensityDPres;
   };
   
   virtual doublereal dGetDensityDTemp(void) const {
      return 0.;
   };
   virtual doublereal dGetDensityDTemp(const doublereal& /* dPres */ ) const {
      return 0.;
   };
   virtual doublereal dGetDensityDTemp(const doublereal& /* dPres */ , const doublereal& /* dTemp */ ) const {
      return 0.;
   };
   
   doublereal dGetViscosity(void) const {
      return dViscosity;
   };
   
   doublereal dGetViscosity(const doublereal& /* dPres */ ) const {
      return dViscosity;
   }; 
   
   doublereal dGetViscosity(const doublereal& /* dPres */ , const doublereal& /* dTemp */ ) const {
      return dViscosity;
   };   
};

/* LinearCompressibleHydraulicFluid - end */


/* LinearCompressibleTHydraulicFluid - begin */

class LinearCompressibleTHydraulicFluid 
: public LinearCompressibleHydraulicFluid {
 protected:
   const doublereal dAlpha; /* mettere il nome vero */
   const doublereal dDensityDTemp;
   
 public:
   LinearCompressibleTHydraulicFluid(unsigned int Label, 
				     const doublereal& dDensity,
				     const doublereal& dBeta,
				     const doublereal& dPres0,
				     const doublereal& dAlpha,
				     const doublereal& dTemp0,
				     const doublereal& dViscosity = 0.)
     : LinearCompressibleHydraulicFluid(Label, dDensity, dBeta, dPres0, dViscosity, dTemp0),
   dAlpha(dAlpha),
   dDensityDTemp(0.) {
      (doublereal&)dDensityDTemp = -dDensity*dAlpha;
   };
   
   LinearCompressibleTHydraulicFluid(const LinearCompressibleTHydraulicFluid& HF)
     : LinearCompressibleHydraulicFluid(HF),
   dAlpha(HF.dAlpha),
      dDensityDTemp(0.) {
      (doublereal&)dDensityDTemp = -dDensity*dAlpha;
   };
   
   ~LinearCompressibleTHydraulicFluid() {
      NO_OP;
   };
   
   /* crea una copia dell'HF */
   HydraulicFluid* pCopy(void) const {
      HydraulicFluid* pHF = NULL;
      SAFENEWWITHCONSTRUCTOR(pHF, 
			     LinearCompressibleTHydraulicFluid,
			     LinearCompressibleTHydraulicFluid(*this));
      return pHF;
   };
   
   /* densita' */
   doublereal dGetDensity(const doublereal& dPres, const doublereal& dTemp) const {
      return LinearCompressibleHydraulicFluid::dGetDensity(dPres)+(dTemp-dTemp0)*dDensityDTemp;
   };
   

   doublereal dGetDensityDTemp(void) const {
      return dDensityDTemp;
   };
   doublereal dGetDensityDTemp(const doublereal& /* dPres */ ) const {
      return dDensityDTemp;
   };
   doublereal dGetDensityDTemp(const doublereal& /* dPres */ , const doublereal& /* dTemp */ ) const {
      return dDensityDTemp;
   };
};

/* LinearCompressibleTHydraulicFluid - end */


/* SuperHydraulicFluid - begin */

static const doublereal a = 2.e-5;
static const doublereal dPresRif = .2*101325.;
static const doublereal dRhoZero = 977.*1.e-2;

class SuperHydraulicFluid : public HydraulicFluid {
 protected:
   const doublereal dDensity;
   const doublereal dBeta;
   const doublereal dDensityDPres;
   const doublereal dViscosity;
   
 public:
   SuperHydraulicFluid(unsigned int Label, 
		       const doublereal& dDens,
		       const doublereal& dB,
		       const doublereal& dPr0,
		       const doublereal& dVisc = 0.,
		       const doublereal dTp0 = -1.)
     : HydraulicFluid(Label, dPr0, dTp0), 
   dDensity(dDens), 
   dBeta(dB),
   dDensityDPres(0.),
   dViscosity(dVisc) {
      if (dB == 0.) {
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
      (doublereal&)dDensityDPres = dDens/dB;
   };
   
   SuperHydraulicFluid(const SuperHydraulicFluid& HF)
     : HydraulicFluid(HF), 
   dDensity(HF.dDensity), 
   dBeta(HF.dBeta),
   dDensityDPres(0.),
   dViscosity(HF.dViscosity) {
      if (dBeta == 0.) {
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
      (doublereal&)dDensityDPres = dDensity/dBeta;
   };
   
   virtual ~SuperHydraulicFluid() {
      NO_OP;
   };
   
   /* crea una copia dell'HF */
   HydraulicFluid* pCopy(void) const {
      HydraulicFluid* pHF = NULL;
      SAFENEWWITHCONSTRUCTOR(pHF,
			     SuperHydraulicFluid,
			     SuperHydraulicFluid(*this));
      return pHF;
   };

   enum { DEFAULT, ZERO, SMALL, INFTY, TANH };
   static const int type = TANH;
   
   /* densita' */
   virtual doublereal dGetDensity(void) const {
      return dDensity;
   };
   
   virtual doublereal dGetDensity(const doublereal& dPres) const {     
      switch (type) {	 
       case DEFAULT:
       default:
	 return dDensity+(dPres-dPres0)*dDensityDPres;
       case ZERO:
	 if (dPres <= dPres0) {
	    return dDensity;
	 }      
	 return dDensity+(dPres-dPres0)*dDensityDPres;
	 
       case SMALL:
	 if (dPres <= 0.) {
	    return dDensity-dPres0*dDensityDPres+dPres*dDensityDPres*1.e-3;
	 }      
	 return dDensity+(dPres-dPres0)*dDensityDPres;
	 
       case INFTY:
	 if (dPres <= 0.) {
	    return dDensity-dPres0*dDensityDPres+dPres*1.e6;
	 }      
	 return dDensity+(dPres-dPres0)*dDensityDPres;
       case TANH:
	 doublereal d = ::dRhoZero+dDensity*.5*(1. + tanh(::a*(dPres-::dPresRif)));
	 if (dPres > ::dPresRif) {
	    d += (dPres-::dPresRif)*dDensityDPres;
	 }
	 return d;
      }
   };
   virtual doublereal dGetDensity(const doublereal& dPres,
				  const doublereal& /* dTemp */ ) const {
      return dGetDensity(dPres);
   };
   

   virtual doublereal dGetDensityDPres(void) const {
      return dDensityDPres;
   };
   
   virtual doublereal dGetDensityDPres(const doublereal& dPres) const {
      switch (type) {
       case DEFAULT: 
       default:
	 return dDensityDPres;
       case ZERO:
	 if (dPres <= 0.) {	    
	    return 0.;
	 }
	 return dDensityDPres;
       case SMALL:
	 if (dPres <= 0.) {	    
	    return dDensityDPres*1.e-3;
	 }
	 return dDensityDPres;
       case INFTY:
	 if (dPres <= 0.) {	   
	    return 1.e6;
	 }
	 return dDensityDPres;
       case TANH:
	 doublereal d = ::a*dDensity*.5/pow(cosh(::a*(dPres-::dPresRif)), 2);
	 if (dPres > ::dPresRif) {
	    d += dDensityDPres;
	 }
	 return d;
      }
   };
   virtual doublereal dGetDensityDPres(const doublereal& dPres, const doublereal& /* dTemp */ ) const {
      return dGetDensityDPres(dPres);
   };
   
   virtual doublereal dGetDensityDTemp(void) const {
      return 0.;
   };
   virtual doublereal dGetDensityDTemp(const doublereal& /* dPres */ ) const {
      return 0.;
   };
   virtual doublereal dGetDensityDTemp(const doublereal& /* dPres */ , const doublereal& /* dTemp */ ) const {
      return 0.;
   };
   
   doublereal dGetViscosity(void) const {
      return dViscosity;
   };
   
   doublereal dGetViscosity(const doublereal& dPres) const {
      if (dPres < ::dPresRif) {
	 return 1.e-6*dViscosity;
      }
      return dViscosity;
   }; 
   
   doublereal dGetViscosity(const doublereal& dPres, const doublereal& /* dTemp */ ) const {
      return dGetViscosity(dPres);
   };   
};

/* SuperHydraulicFluid - end */

/* ExpHydraulicFluid - begin */

class ExpHydraulicFluid : public HydraulicFluid {
 protected:
   const doublereal dDensity;
   const doublereal dBeta;
   const doublereal dViscosity;
   const doublereal dPsat;
   
 public:
   ExpHydraulicFluid(unsigned int Label, 
			const doublereal& dDens,
			const doublereal& dB,
			const doublereal& dPr0,
			const doublereal& dPs,
			const doublereal& dVisc = 0.,
			const doublereal dTp0 = -1.)
     : HydraulicFluid(Label, dPr0, dTp0), 
   dDensity(dDens), 
   dBeta(dB),
   dViscosity(dVisc),
   dPsat(dPs) {
      if (dB == 0.) {
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
   };
   
   ExpHydraulicFluid(const ExpHydraulicFluid& HF)
     : HydraulicFluid(HF), 
   dDensity(HF.dDensity), 
   dBeta(HF.dBeta),
   dViscosity(HF.dViscosity),
   dPsat(HF.dPsat) {
      if (dBeta == 0.) {
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }
   };
   
   virtual ~ExpHydraulicFluid() {
      NO_OP;
   };
   
   /* crea una copia dell'HF */
   HydraulicFluid* pCopy(void) const {
      HydraulicFluid* pHF = NULL;
      SAFENEWWITHCONSTRUCTOR(pHF,
			     ExpHydraulicFluid,
			     ExpHydraulicFluid(*this));
      return pHF;
   };

   /* densita' */
   virtual doublereal dGetDensity(void) const {
      return dDensity;
   };
   
   virtual doublereal dGetDensity(const doublereal& dPres) const {
      /*
       * Inspired by AMESim's simple saturating fluid 
       */
      if (dPres >= dPsat) {
	 return dDensity*exp((dPres-dPres0)/dBeta);
      } else {
	 return dDensity*exp(1000.*(dPres-dPres0)/dBeta);
      }
   };
   virtual doublereal dGetDensity(const doublereal& dPres,
				  const doublereal& /* dTemp */ ) const {
      return dGetDensity(dPres);
   };
   

   virtual doublereal dGetDensityDPres(void) const {
      return dDensity/dBeta;
   };
   
   virtual doublereal dGetDensityDPres(const doublereal& dPres) const {
      if (dPres >= dPsat) {
	 return dDensity/dBeta*exp((dPres-dPres0)/dBeta);
      } else {
	 return dDensity*1000./dBeta*exp(1000.*(dPres-dPres0)/dBeta);
      }
   };
   
   virtual doublereal dGetDensityDPres(const doublereal& dPres, const doublereal& /* dTemp */ ) const {
      return dGetDensityDPres(dPres);
   };
   
   virtual doublereal dGetDensityDTemp(void) const {
      return 0.;
   };
   virtual doublereal dGetDensityDTemp(const doublereal& /* dPres */ ) const {
      return 0.;
   };
   virtual doublereal dGetDensityDTemp(const doublereal& /* dPres */ , const doublereal& /* dTemp */ ) const {
      return 0.;
   };
   
   doublereal dGetViscosity(void) const {
      return dViscosity;
   };
   
   doublereal dGetViscosity(const doublereal& dPres) const {
      if (dPres >= dPsat) {
	 return dViscosity;
      } else {
	 return 0.;
      }
   }; 
   
   doublereal dGetViscosity(const doublereal& dPres, const doublereal& /* dTemp */ ) const {
      return dGetViscosity(dPres);
   };
};

/* ExpHydraulicFluid - end */

#endif /* HFLUID__H */
