/*

MBDyn (C) is a multibody analysis code. 
http://www.mbdyn.org

Copyright (C) 1996-2007

Pierangelo Masarati	<masarati@aero.polimi.it>
Paolo Mantegazza	<mantegazza@aero.polimi.it>

Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
via La Masa, 34 - 20156 Milano, Italy
http://www.aero.polimi.it

Changing this copyright notice is forbidden.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


------------------------------------------------------------------------------

ADAMS2MBDyn (C) is a translator from ADAMS/View models in adm format
into raw MBDyn input files.

Copyright (C) 1999-2007
Leonardo Cassan		<lcassan@tiscalinet.it>

*/


#ifndef IC_H
#define IC_H

#include <common.h>
#include <output.h>

struct s_ic : public s_card {
   double Aerror;
   Angle Alimit;
   int Amaxit;
   Boolean* Apattern;
   double Error;
   int Maxit;
   Boolean* Pattern;
   double Tlimit;
   double Verror;
   int ncoord_pattern;
   int ncoord_apattern;
   
   Boolean _Aerror, _Alimit, _Amaxit, _Apattern,
           _Error, _Maxit, _Pattern, _Tlimit, _Verror;
   
   s_ic(void);
   ~s_ic(void);
   
   ostream& Print (ostream& out) const;
   void Translate (ostream& out);
   inline const char* const Gettype(void) const;
   Boolean Test();
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch (param) {
	 case AERROR : Set_Param (Aerror,value,_Aerror); break;
	 case ALIMIT: Set_Param(Alimit,value,_Alimit); break;
	 case AMAXIT: Set_Param(Amaxit,value,_Amaxit); break;
	 case APATTERN: Set_Array (Apattern,value,ncoord_apattern,_Apattern); break;
	 case ERROR: Set_Param(Error,value,_Error); break;
	 case MAXIT: Set_Param(Maxit,value,_Maxit); break;
	 case PATTERN: Set_Array (Pattern,value,ncoord_pattern,_Pattern); break;
	 case TLIMIT: Set_Param(Tlimit,value,_Tlimit); break;
	 case VERROR: Set_Param(Verror,value,_Verror); break;
	 default: out_error (3,Find_Token(param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
	return;
     }
};

#endif
