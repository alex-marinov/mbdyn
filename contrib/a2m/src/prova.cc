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




#include <matrix.h>
#include <iostream>
#include <math.h>

Vec3 a(1,0,0),b(50,50,0);
double k1,k2;
Eye PR(3);
Matrix B (3,1);

main ()
{
   const char indent[]="       ";
   cout << "Hello!" << endl;
   cout << "Vectors are:" << a << ", " << b << endl;
   k1=a*b;
   cout << "Prodotto scalare:" << k1 << endl;
   cout << "Coseno angolo:" << Cos_b2Vector(a,b) << endl;
   cout << "Angolo:" << (57.3*acos(Cos_b2Vector(a,b))) << endl;
   Vector Q(3);
   Q=a^b;
   cout << "Prodotto vettore:" << Q << endl;
   cout << "Risultato dell'interpolazione" << endl;
   cout << a << endl << b << endl;
   cout << Interp(a,b) << endl;
   cout << (Interp(a,b)/25.5) << endl;
   cout << "",PR.Write (cout,",","\n",indent);
   cout << endl;
}
