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

/* Genels */

#ifndef GENEL_H
#define GENEL_H

#include "elem.h"
#include "node.h"

extern const char* psGenelNames[];


/* Genel - begin */

class Genel : virtual public Elem, public ElemWithDofs {
 public:
   /* Tipi di Genel */
   enum Type {
      UNKNOWN = -1,
	SWASHPLATE = 0,
	ROTORTRIM,
	CLAMP,
	DISTANCE,
	SPRING,
	SPRINGSUPPORT,
	CROSSSPRINGSUPPORT,
	SPRINGDAMPER,
	SPRINGDAMPERSUPPORT,
	CROSSSPRINGDAMPERSUPPORT,
	MASS,
	SCALARFILTER,
	STATESPACESISO,
	STATESPACEMIMO,
	
	LASTGENELTYPE
   };

 public:
   Genel(unsigned int uL, const DofOwner* pDO, flag fOut);
   virtual ~Genel(void);
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   /* Tipo dell'elemento (usato per debug ecc.) */
   virtual Elem::Type GetElemType(void) const {
      return Elem::GENEL; 
   };
   
   /* Tipo di Genel */
   virtual Genel::Type GetGenelType(void) const = 0;
};
   
/* Genel - end */

   
class DataManager;
class MBDynParser;

extern Elem* ReadGenel(DataManager* pDM, 
		       MBDynParser& HP,
		       const DofOwner* pDO, 
		       unsigned int uLabel);

#endif /* GENEL_H */

