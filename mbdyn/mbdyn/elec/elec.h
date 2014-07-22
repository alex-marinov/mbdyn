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

/* elementi elettrici, tipo: Elem::ELECTRIC */

#ifndef ELEC_H
#define ELEC_H

#include "elem.h"
      
extern const char* psElectricNames[];

/* Electric - begin */

class Electric : virtual public Elem, public ElemWithDofs {
 public:
   /* Tipi di elementi elettrici */
   enum Type {
      UNKNOWN = -1,
      
      ACCELEROMETER = 0,
      DISPLACEMENT,
      DISCRETECONTROL,
      MOTOR,
      
      LASTELECTRICTYPE
   };

 public:
   Electric(unsigned int uL,
	    const DofOwner* pDO, flag fOut);
   virtual ~Electric(void);
   
   /* Contributo al file di restart 
    * (Nota: e' incompleta, deve essere chiamata dalla funzione corrispndente
    * relativa alla classe derivata */
   virtual std::ostream& Restart(std::ostream& out) const;

   /* Tipo dell'elemento (usato solo per debug ecc.) */
   virtual Elem::Type GetElemType(void) const;

   /* Tipo di elemento elettrico (usato solo per debug ecc.) */
   virtual Electric::Type GetElectricType(void) const = 0;
};

/* Electric - end */

class DataManager;
class MBDynParser;

extern Elem* ReadElectric(DataManager* pDM,
			  MBDynParser& HP, 
			  const DofOwner* pDO, 
			  unsigned int uLabel);

#endif /* ELEC_H */

