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

#ifndef OUTPUT_H
#define OUTPUT_H

// FUNZIONI NECESSARIE ALLA CORRETTA VISUALIZZAZIONE DELL'OUTPUT

#include <defs.h>
#include <iostream>

ostream& operator << (ostream& out, const Joint&);
ostream& operator << (ostream& out, const Friction&) ;
ostream& operator << (ostream& out, const Joint_Primitive&);
ostream& operator << (ostream& out, const Direction_Mode&);
ostream& operator << (ostream& out, const Boolean&);
ostream& operator << (ostream& out, const coord_type&);

template <class T>
ostream& outvec (ostream& out, T* data, int idx)
{
   for (int k=0; k<idx; k++)
     out << data[k] << " ";
   return out;
}

#endif
