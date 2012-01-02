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


#ifndef PRESELEM_H
#define PRESELEM_H

#include "elem.h"
#include "presnode.h"
#include "drive.h"    /* per parametri variabili */
#include "strnode.h"  /* per attuatore */

#include "hfluid.h"


extern const char* psHydraulicNames[];


/* HydraulicElem - begin */

class HydraulicElem : virtual public Elem, public ElemWithDofs {
 public:
   /* Tipi di elementi idraulici */
   enum Type {
      UNKNOWN = -1,
	
	MINOR_LOSS = 0,
	THREEWAYMINORLOSS,
	CONTROL_VALVE,
	DYNAMIC_CONTROL_VALVE,
        PRESSURE_FLOW_CONTROL_VALVE,
        PRESSURE_VALVE,
	FLOW_VALVE,
	ORIFICE,
	ACCUMULATOR,
        TANK,
        PIPE,
        DYNAMIC_PIPE,
	HYDRAULIC_ACTUATOR,
        ACTUATOR,
      
	LASTHYDRAULICTYPE
   };
 
 protected:
   HydraulicFluid* HF;
   
 public:
   HydraulicElem(unsigned int uL, const DofOwner* pDO, 
		 HydraulicFluid* hf, flag fOut);
   virtual ~HydraulicElem(void);
   
   /* Tipo dell'elemento (usato per debug ecc.) */
   virtual Elem::Type GetElemType(void) const;
   
   /* Contributo al file di restart 
    * (Nota: e' incompleta, deve essere chiamata dalla funzione corrispndente
    * relativa alla classe derivata */
   virtual std::ostream& Restart(std::ostream& out) const;

   /* Tipo di elemento elettrico (usato solo per debug ecc.) */
   virtual HydraulicElem::Type GetHydraulicType(void) const = 0;
};

/* HydraulicElem - end */


class DataManager;
class MBDynParser;

extern Elem* ReadHydraulicElem(DataManager* pDM,
			       MBDynParser& HP, 
			       const DofOwner* pDO, 
			       unsigned int uLabel);


#endif /* PRESELEM_H */
