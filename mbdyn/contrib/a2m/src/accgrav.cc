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


//                              ACCGRAV.CC                                    

#include <accgrav.h>

extern MBDyn_deck MBAccs;

s_accgrav::s_accgrav (void)   : Igrav(0),Jgrav(0),Kgrav(0),
                               _Igrav(N),_Jgrav(N),_Kgrav(N)
                               {}

Boolean s_accgrav::Test()
{
   /* non necessita di controlli incrociati */
   return N;
}

inline const char* const s_accgrav::Gettype (void) const
{
   return "ACCGRAV";
}

ostream& s_accgrav::Print (ostream& out) const
{
   out << endl;
   out << "ACCGRAV" << endl
     << "     I COMPONENT [" << _Igrav << "] = " << Igrav << endl
     << "     J COMPONENT [" << _Jgrav << "] = " << Jgrav << endl
     << "     K COMPONENT [" << _Kgrav << "] = " << Kgrav << endl;
   return out;
}

void s_accgrav::Translate(ostream& out)
{
   Vec3 G(Igrav,Jgrav,Kgrav);
   char* comment=new char [255];
   /* La gravità è espressa nel sistema globale.
    * Il template drive_caller di tipo Vec3 è il risultato
    * di una moltiplicazione del vettore Accellerazione per
    * la costante 1 */

   /* crea un drive costante con coefficiente unitario */
   MBDyn_drive_CONST* Accdrive = new MBDyn_drive_CONST (1);
   sprintf (comment,"/* drive */ ");
   Accdrive->Remark(comment);
   
   MBDyn_gravity* GRAVITY = new MBDyn_gravity (G,Accdrive);
   sprintf (comment,"Gravity accelleration in MBDyn format [Adams ACCGRAV]");
   
   GRAVITY->Remark(comment);
   
   MBAccs.insert (MBDyn_entry(999,(MBDyn_card*) GRAVITY));
   return;
}
