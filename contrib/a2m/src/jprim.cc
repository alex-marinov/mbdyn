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

//                                 JPRIM.CC                                   

#include <jprim.h>
#include <output.h>

s_jprim::s_jprim (void)       : Node1(0),Node2(0),Type(_INPLANE),
                                _Node1(N),_Node2(N),_Type(N)
                                {}

s_jprim::~s_jprim (void)      {}


Boolean s_jprim::Test()
{
   const int err_before=nerr;
   // VERIFICA DELL'ESISTENZA DI ENTRAMBI GLI ID DEI MARKER I E J
   if (_Node1==N) out_error (27,"ID I");
   if (_Node2==N) out_error (27,"ID J");
   // VERIFICA CHE SIA STATO SPECIFICATO IL TIPO DI PRIMITIVA
   if (_Type==N) out_error (28,"");
   if (err_before != nerr) return Y; else return N;
}

inline const char* const s_jprim::Gettype (void) const
{
   return "JPRIM";
}

ostream& s_jprim::Print (ostream& out) const
{
   out << endl;
   out << "JPRIM:" << label << "      TYPE:" << Type << endl;
   out << "     Node1 [" << _Node1 << "] = " << Node1 << endl
       << "     Node2 [" << _Node2 << "] = " << Node2 << endl;
   return out;
}

void s_jprim::Translate (ostream& out)
{
   return;
}


