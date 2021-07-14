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

//                                  IC.CC                                      

#include <ic.h>

s_ic::s_ic (void)  :   Aerror(1.0E-4),Alimit(30,DEGREE),Amaxit(25),
                       Apattern(new Boolean[10] ), Pattern(new Boolean[10]),
                       ncoord_pattern(0),ncoord_apattern(0),
                       Error(1.0E-10),Maxit(25),Tlimit(1.0E10),Verror(1E-4),
                       _Aerror(N),_Alimit(N),_Amaxit(N),_Error(N),_Maxit(N),
                       _Tlimit(N),_Verror(N),_Apattern(N),_Pattern(N)
{
   for (int i=0; i<10; i++) {
      Apattern[i]=Y;
      Pattern[i]=Y;
   }
}

s_ic::~s_ic(void)      {}

inline const char* const s_ic::Gettype (void) const
{
   return "IC";
}

Boolean s_ic::Test()
{
   /* Non necessita di controlli incrociati */
   return N;
}

ostream& s_ic::Print (ostream& out) const
{
   out << endl;
   out << "IC:" << endl;
   out << "     Aerror [" << _Aerror << "] = " << Aerror << endl
       << "     Alimit [" << _Alimit << "] = " << Alimit << endl
       << "     Amaxit [" << _Amaxit << "] = " << Amaxit << endl
       << "     Error  [" << _Error << "] = " << Error << endl
       << "     Maxit  [" << _Maxit << "] = " << Maxit << endl
       << "     Tlimit [" << _Tlimit << "] = " << Tlimit << endl
       << "     Verror [" << _Verror << "] = " << Verror << endl
       << "     Apattern [" << _Apattern << "] = ";
   if (_Apattern==Y) outvec (out,Apattern,ncoord_apattern);
   out << endl;
   out << "     Pattern [" << _Pattern << "] = ";
   if (_Pattern==Y) outvec(out,Pattern,ncoord_pattern);
   out << endl;
   out << endl;
   return out;
}

void s_ic::Translate (ostream& out)
{
   return;
}



