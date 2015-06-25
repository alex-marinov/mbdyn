/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include "bulk.h"
#include "dataman.h"

/* Legge un elemento bulk */
   
Elem* ReadBulk(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
   DEBUGCOUTFNAME("ReadBulk");
   
   const char* sKeyWords[] = {
      "spring" "support",
      NULL
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
	SPRINGSUPPORT = 0,	
	LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K(HP, sKeyWords);
   
   /* lettura del tipo di elemento elettrico */   
   KeyWords CurrKeyWord = KeyWords(HP.GetWord());
   
#ifdef DEBUG   
   if (CurrKeyWord >= 0) {      
      std::cout << "bulk element type: " 
	<< sKeyWords[CurrKeyWord] << std::endl;
   }   
#endif   

   Elem* pEl = NULL;
   
   switch (CurrKeyWord) {
      /*  */
    case SPRINGSUPPORT: {       
       ScalarDof SD = ReadScalarDof(pDM, HP, true, false);
       DriveCaller* pDC = HP.GetDriveCaller();
       flag fOut = pDM->fReadOutput(HP, Elem::BULK);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      BulkSpringSupport,
			      BulkSpringSupport(uLabel, pDC, SD, fOut));
       
       break;
    }
      
      /* Aggiungere altri elementi elettrici */
      
    default: {
       silent_cerr("unknown bulk element type in bulk element " << uLabel
	 << " at line " << HP.GetLineData() << std::endl);
       throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
    }	
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.IsArg()) {
      silent_cerr("semicolon expected at line " 
	      << HP.GetLineData() << std::endl);
      throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
   }   
   
   return pEl;
} /* End of ReadBulk() */
