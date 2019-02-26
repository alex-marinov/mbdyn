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


//                                EQUILIBRIUM                                 

// Definizioni per la card Equilibrium

#include <equilibrium.h>

s_equilibrium::s_equilibrium (void) : Alimit(0,DEGREE),Error(1.0E-4),
                                      Imbalance(1.0E-4),Maxit(25),
                                      Stability(1.0E-5),Tlimit(1.0E10),
                                      Pattern(new Boolean[10]),
                                      ncoord_pattern(0)
{
   for (int i=0; i<10; i++) Pattern[i]=Y;
}

s_equilibrium::~s_equilibrium (void)  {
}

Boolean s_equilibrium::Test()
{
   /* Non si rende necessario l'utilizzo di controlli incrociati */
   return N;
}

inline const char* const s_equilibrium::Gettype(void) const
{
   return "EQUILIBRIUM";
}

ostream& s_equilibrium::Print (ostream& out) const
{
   out << endl;
   out << "EQUILIBRIUM:" << endl;
   out << "     " << "Alimit [" << _Alimit << "] = " << Alimit << endl;
   out << "     " << "Error [" << _Error << "] = " << Error << endl
       << "     " << "Imbalance [" << _Imbalance << "] = "
       << Imbalance << endl 
       << "     " << "Max iterations [" << _Maxit << "] = "
       << Maxit << endl
       << "     " << "Stability [" << _Stability << "] = "
       << Stability << endl
       << "     " << "Tlimit [" << _Tlimit << "] = " 
       << Tlimit << endl;
   out << "     " << "Pattern [" << _Pattern << "] = ";
   if (_Pattern==Y) outvec(out,Pattern,ncoord_pattern);
   out << endl << endl;
   return out;
}

void s_equilibrium::Translate (ostream& out)
{
   return;
}
