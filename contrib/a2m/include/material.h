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

#ifndef MATERIAL_H
#define MATERIAL_H

#include <common.h>

struct s_material : public s_card {
   double Young;
   double Density;
   double P_ratio;
   char *Name;
   
   Boolean _Young,_Density,_Name,_P_ratio;
   
   s_material(void);
   ~s_material (void);
   ostream& Print (ostream& out) const;
   void Translate(ostream& out);
   inline const char* const Gettype (void) const;
   Boolean Test();
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch (param) {
	 case YOUNGS_MODULUS : Set_Param (Young,value,_Young); break;
	 case DENSITY : Set_Param (Density, value, _Density); break;
	 case NAME : Set_Param (Name, value, _Name); break;
	 case POISSONS_RATIO : Set_Param (P_ratio, value, _P_ratio); break; 
	 default : out_error (3,Find_Token(param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
     }
};

#endif
