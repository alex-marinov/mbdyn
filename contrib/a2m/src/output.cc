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

//                                OUTPUT.CC                                   


#include <output.h>


ostream& operator << (ostream& out, const coord_type& R)
{
   switch(R) {
    case _XC: out << "X"; break;
    case _YC: out << "Y"; break;
    case _ZC: out << "Z"; break;
    case _PSIC : out << "PSI"; break;
    case _THETAC : out << "THETA"; break;
    case _PHIC : out << "PHI"; break;
    default: out << "NULL"; break;
   }
   return out;
}

ostream& operator << (ostream& out, const Joint& R)
{
   switch (R) {
    case _CYLINDRICAL   : out << "CYLINDRICAL"; break;
    case _CONVEL	: out << "CONVEL"; break;
    case _REVOLUTE	: out << "REVOLUTE"; break;
    case _SPHERICAL	: out << "SPHERICAL"; break;
    case _SCREW         : out << "SCREW"; break;
    case _FIXED         : out << "FIXED"; break;
    case _HOOKE         : out << "HOOKE"; break;
    case _TRANSLATIONAL : out << "TRANSLATIONAL"; break;
    case _UNIVERSAL     : out << "UNIVERSAL"; break;
    case _RACKPIN       : out << "RACKPIN"; break;
    case _PLANAR        : out << "PLANAR"; break;
    default		: out << "UNCLASSIFIED."; break;
   }
   return out;
}

ostream& operator << (ostream& out, const Boolean& R)
{
   if (R==Y) out << "TRUE"; else out << "FALSE";
   return out;
}

ostream& operator << (ostream& out, const Friction& R)
{
   switch(R) {
    case _OFF         : out << "OFF"; break;
    case _ON          : out << "ON"; break;
    case _PRELOAD_ONLY: out << "PRELOAD_ONLY"; break;
   }
   return out;
}

ostream& operator << (ostream& out, const Joint_Primitive& R)
{
   switch (R) {
    case _ATPOINT: out << "ATPOINT"; break;
    case _INLINE: out << "INLINE"; break;
    case _INPLANE: out << "INPLANE"; break;
    case _ORIENTATION: out << "ORIENTATION"; break;
    case _PARALLEL_AXES: out << "PARALLEL AXES"; break;
    case _PERPENDICULAR: out << "PERPENDICULAR"; break;
   }
   return out;
}

ostream& operator <<  (ostream& out, const Direction_Mode& R)
{
   switch (R) {
    case _ROTATION: out << "ROTATION"; break;
    case _TRANSLATION: out << "TRANSLATION"; break;
   }
   return out;
}


