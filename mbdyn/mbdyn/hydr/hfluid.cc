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

/* 
 * Copyright 1999-2000 Lamberto Puggelli <puggelli@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

#include <mbconfig.h>

#include <dataman.h>
#include <hfluid.h>
#include <hfluid_.h>

HydraulicFluid*
ReadHydraulicFluid(MBDynParser& HP,
		   unsigned int uLabel)
{
   /* tipi di fluidi */
   const char* sKeyWords[] = {
      "uncompressible",
      "linear" "compressible",
      "linear" "thermal" "compressible",
      "super",
      "amesim",
      
      NULL
   };
   
   enum KeyWords {
      UNKNOWN = -1,
      
      UNCOMPRESSIBLE = 0,
      LINEARCOMPRESSIBLE,
      LINEARTHERMALCOMPRESSIBLE,
      SUPER,
      AMESIM,
      
      LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* attacca la tabella al parser */
   HP.PutKeyTable(K);
   
   /* lettura del tipo di drive */   
   KeyWords CurrKeyWord;
   if ((CurrKeyWord = KeyWords(HP.IsKeyWord())) == UNKNOWN) {
      CurrKeyWord = UNCOMPRESSIBLE;
   }
   
   HydraulicFluid* pHF = NULL;
   
   switch (CurrKeyWord) {
    case UNCOMPRESSIBLE: {
       
       doublereal dDensity(0.);
       if (HP.IsKeyWord("density")) {
	  dDensity = HP.GetReal();
       }
       if (dDensity < DBL_EPSILON) {
	  cerr << "line " << HP.GetLineData()
	    << ": illegal density " << dDensity
	    << " for hydraulic fluid " << uLabel << endl;
       }
       
       doublereal dViscosity(0.);
       if (HP.IsKeyWord("viscosity")) {
	  dViscosity = HP.GetReal();
       }
        if (dViscosity < DBL_EPSILON) {
		  cerr << "line " << HP.GetLineData() 
	            << ": illegal viscosity " << dViscosity
		    << endl;
		  THROW(ErrGeneric());
	       }
       
       doublereal dPres0(-1.);
       if (HP.IsKeyWord("pressure")) {
	  dPres0 = HP.GetReal();
	  if (dPres0 < 0.) {
	     cerr << "line " << HP.GetLineData() 
	       << ": illegal reference pressure " << dPres0 
	       << " for hydraulic fluid " << uLabel << endl;
	  }	       
       }
       
       doublereal dTemp0(-1.);
       if (HP.IsKeyWord("temperature")) {
	  dTemp0 = HP.GetReal();
	  if (dTemp0 < 0.) {
	     cerr << "line " << HP.GetLineData()
	       << ": illegal reference temperature " << dTemp0
	       << " for hydraulic fluid " << uLabel << endl;
	  }	       
       }       
       
       SAFENEWWITHCONSTRUCTOR(pHF,
			      UncompressibleHydraulicFluid,
			      UncompressibleHydraulicFluid(uLabel, 
							   dDensity, 
							   dViscosity,
							   dPres0,
							   dTemp0),
			      DMmm);
       break;
    }
      
    case LINEARCOMPRESSIBLE: 
    case SUPER:
    case AMESIM: {
       doublereal dDensity(0.);
       doublereal dBeta(0.);
       doublereal dPres0(0.);
       if (HP.IsKeyWord("density")) {
	  dDensity = HP.GetReal();
	   
	  if (HP.IsKeyWord("soundcelerity")) {
	     doublereal sound = HP.GetReal();
	     if (sound < DBL_EPSILON) {
		cerr << "line " << HP.GetLineData() << ": illegal sound velocity " << sound
		  << endl;
		THROW(ErrGeneric());		
	     } 
	     dBeta = sound*sound*dDensity;
	  } else
	    { 
	       dBeta = HP.GetReal();
	       if (fabs(dBeta) < DBL_EPSILON) {
		  cerr << "line " << HP.GetLineData() << ": illegal bulk modulus " << dBeta 
		    << endl;
		  THROW(ErrGeneric());
	       }
	    }
	   
	  if (dDensity < DBL_EPSILON) {
	     cerr << "line " << HP.GetLineData()
	       << ": illegal density " << dDensity
	       << " for hydraulic fluid " << uLabel << endl;
	  }
	  
	  dPres0 = HP.GetReal();
	  if (dPres0 < 0.) {
	     cerr << "line " << HP.GetLineData() 
	       << ": illegal reference pressure " << dPres0 
	       << " for hydraulic fluid " << uLabel << endl;
	  }	       
       }
       
       doublereal dViscosity(0.);
       if (HP.IsKeyWord("viscosity")) {
	  dViscosity = HP.GetReal();
       }
        
       if (dViscosity < DBL_EPSILON) {
	  cerr << "line " << HP.GetLineData() 
	    << ": illegal viscosity " << dViscosity
	    << endl;
	  THROW(ErrGeneric());
       }
       
       doublereal dTemp0(-1.);
       if (HP.IsKeyWord("temperature")) {
	  dTemp0 = HP.GetReal();
	  if (dTemp0 < 0.) {
	     cerr << "line " << HP.GetLineData()
	       << ": illegal reference temperature " << dTemp0
	       << " for hydraulic fluid " << uLabel << endl;
	  }	       
       }
       
       switch (CurrKeyWord) {
	case LINEARCOMPRESSIBLE: 
	  SAFENEWWITHCONSTRUCTOR(pHF,
				 LinearCompressibleHydraulicFluid,
				 LinearCompressibleHydraulicFluid(uLabel, 
								  dDensity, 
								  dBeta,
								  dPres0,
								  dViscosity,
								  dTemp0),
				 DMmm);
	  break;
	case SUPER:
	  SAFENEWWITHCONSTRUCTOR(pHF,
				 SuperHydraulicFluid,
				 SuperHydraulicFluid(uLabel, 
						     dDensity, 
						     dBeta,
						     dPres0,
						     dViscosity,
						     dTemp0),
				 DMmm);
	  break;
	case AMESIM:
	  doublereal dPsat = HP.GetReal();
	  
	  SAFENEWWITHCONSTRUCTOR(pHF,
				 AMESimHydraulicFluid,
				 AMESimHydraulicFluid(uLabel, 
						      dDensity, 
						      dBeta,
						      dPres0,
						      dPsat,
						      dViscosity,
						      dTemp0),
				 DMmm);
	  break;
       }
       
       break;
    }

    case LINEARTHERMALCOMPRESSIBLE: {
       doublereal dDensity(0.);
       doublereal dBeta(0.);
       doublereal dPres0(0.);
       doublereal dAlpha(0.);
       doublereal dTemp0(0.);
       if (HP.IsKeyWord("density")) {
	  dDensity = HP.GetReal();
	 
	  if (HP.IsKeyWord("soundcelerity")) {
	     doublereal sound = HP.GetReal();
	     if (sound < DBL_EPSILON) {
		cerr << "line " << HP.GetLineData() << ": illegal sound velocity " << sound
		  << endl;
		THROW(ErrGeneric());		
	     } 
	     dBeta = sound*sound*dDensity;
	  } else {
	     dBeta = HP.GetReal();
	     if (fabs(dBeta) < DBL_EPSILON) {
		cerr << "line " << HP.GetLineData() << ": illegal bulk modulus " << dBeta << endl;
		THROW(ErrGeneric());
	     }
	  }
	  
	  if (dDensity < DBL_EPSILON) {
	     cerr << "line " << HP.GetLineData()
	       << ": illegal density " << dDensity
	       << " for hydraulic fluid " << uLabel << endl;
	  }
	  
	  dPres0 = HP.GetReal();
	  if (dPres0 < 0.) {
	     cerr << "line " << HP.GetLineData() 
	       << ": illegal reference pressure " << dPres0 
	       << " for hydraulic fluid " << uLabel << endl;
	  }	       

	  dAlpha = HP.GetReal();
	  
	  dTemp0 = HP.GetReal();
	  if (dTemp0 < 0.) {
	     cerr << "line " << HP.GetLineData()
	       << ": illegal reference temperature " << dTemp0
	       << " for hydraulic fluid " << uLabel << endl;
	  }	       
       }
       
       doublereal dViscosity(0.);
       if (HP.IsKeyWord("viscosity")) {
	  dViscosity = HP.GetReal();
       }
       
       if (dViscosity < DBL_EPSILON) {
	  cerr << "line " << HP.GetLineData() 
	    << ": illegal viscosity " << dViscosity
	    << endl;
	  THROW(ErrGeneric());
       }
       
       SAFENEWWITHCONSTRUCTOR(pHF,
			      LinearCompressibleTHydraulicFluid,
			      LinearCompressibleTHydraulicFluid(uLabel,
								dDensity,
								dBeta,
								dPres0,
								dAlpha,
								dTemp0,
								dViscosity),
			      DMmm);
       break;
    }

    default:
      cerr << "line " << HP.GetLineData() << ": unknown hydraulic fluid type" << endl;
      THROW(ErrGeneric());
   }
   
   return pHF;
}
