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

//                                 VARIABLE.CC                                  


#include <variable.h>
#include <output.h>


s_variable::s_variable(void) : Ic(0), Expression (new char[256]),
                               Value(0),
                               _Ic(N),_Expression(N),_Value(N) { }

s_variable::~s_variable(void) {}

inline const char* const s_variable::Gettype(void) const
{
   return "VARIABLE";
}

Boolean s_variable::Test(void)
{
   const int err_before=nerr;
   /* here is the test code */
   if (err_before != nerr) return Y; else return N;
}

ostream& s_variable::Print (ostream& out) const
{
   out << endl;
   out << "VARIABLE: " << label << endl;
   out << "      " << "Ic [" << _Ic << "] = " << Ic << endl
       << "      " << "Expression [" << _Expression << "] = " << Expression
       << endl
       << "      " << "Value [" << _Value << "] = " << Value
       << endl;
   out << endl;
   return out;
}

void s_variable::Translate (ostream& out)
{
   /* Here is the translation code */
   return;
}
