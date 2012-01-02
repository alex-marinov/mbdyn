/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#ifndef HFLUID_H
#define HFLUID_H

#include "ac/f2c.h"

#include "withlab.h"

/* HydraulicFluid - begin */

class HydraulicFluid : public WithLabel {
 protected: 
   const doublereal dPres0;
   const doublereal dTemp0;
   
 public:
   HydraulicFluid(unsigned int Label,
		  const doublereal dPres0 = -1.,
		  const doublereal dTemp0 = -1.) 
     : WithLabel(Label), dPres0(dPres0), dTemp0(dTemp0) {
	NO_OP;
     };
   HydraulicFluid(const HydraulicFluid& HF)
     : WithLabel(HF.GetLabel()), dPres0(HF.dPres0), dTemp0(HF.dTemp0) {
	NO_OP;
     };
   
   virtual ~HydraulicFluid() {
      NO_OP;
   };
   
   /* crea una copia dell'HF */
   virtual HydraulicFluid* pCopy(void) const = 0;

   /* funzioni relative alla densita'; possono essere usate senza argomenti
    * (fluido incomprimibile), con un argomento, o con due (effetto di T) */
   virtual doublereal dGetDensity(void) const = 0;
   virtual doublereal dGetDensity(const doublereal& dPres) const = 0;
   virtual doublereal dGetDensity(const doublereal& dPres,
				  const doublereal& dTemp) const = 0;
   virtual doublereal dGetDensityDPres(void) const = 0;
   virtual doublereal dGetDensityDPres(const doublereal& dPres) const = 0;
   virtual doublereal dGetDensityDPres(const doublereal& dPres,
				       const doublereal& dTemp) const = 0;
   virtual doublereal dGetDensityDTemp(void) const = 0;
   virtual doublereal dGetDensityDTemp(const doublereal& dPres) const = 0;
   virtual doublereal dGetDensityDTemp(const doublereal& dPres,
				       const doublereal& dTemp) const = 0;
   virtual doublereal dGetViscosity(void) const = 0;
   virtual doublereal dGetViscosity(const doublereal& dPres) const = 0;
   virtual doublereal dGetViscosity(const doublereal& dPres,
				    const doublereal& dTemp) const = 0;
   enum Re { UPPER, LOWER };
   virtual inline doublereal dGetRe(Re which);
#if 0
   virtual doublereal dGetRe(Re which, const doublereal& dPres) = 0;
   virtual doublereal dGetRe(Re which, const doublereal& dPres,
		   const doublereal& dTemp) = 0;
#endif
   
   virtual doublereal dGetPres0(void) const {
      ASSERT(dPres0 != -1.);
      return dPres0;
   };
   virtual doublereal dGetTemp0(void) const {
      ASSERT(dTemp0 != -1.);
      return dTemp0;
   };
};

inline doublereal 
HydraulicFluid::dGetRe(HydraulicFluid::Re which)
{
	switch (which) {
	case HydraulicFluid::LOWER:
		return 2000.;
		
	case HydraulicFluid::UPPER:
		return 4000.;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* legge un fluido idraulico; non va usata direttamente, viene invece 
 * gestita dal parser specializzato di MBDyn, MBDynParser (mbpar.h) */
class MBDynParser;
extern HydraulicFluid* ReadHydraulicFluid(MBDynParser& HP, unsigned int uLabel);

/* HydraulicFluid - end */

#endif /* HFLUID_H */

