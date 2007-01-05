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

//                              MATERIAL.CC                                   

#include <material.h>

s_material::s_material (void) : Young(0),Density(0),Name(),P_ratio(0),
                               _Young(N),_Density(N),_Name(N),_P_ratio(N) 
                               { Name=new char[80]; }

s_material::~s_material (void)  {}


Boolean s_material::Test()
{
   /* Rimanda in un secondo tempo */
   return N;
}

inline const char* const s_material::Gettype (void) const
{
   return "MATERIAL";
}

ostream& s_material::Print (ostream& out) const
{
   out << endl;
   out << "MATERIAL:" << label << "       NAME:" << Name << endl
     << "     YOUNG MODULUS [" << _Young << "] = " << Young << endl
     << "     DENSITY [" << _Density << "] = " << Density << endl
     << "     POISSONS RATIO [" << _P_ratio << "] = " << P_ratio << endl;
     return out;
}

void s_material::Translate(ostream& out)
{
   return;
}


