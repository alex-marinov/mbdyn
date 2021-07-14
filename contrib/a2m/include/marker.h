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

#ifndef MARKER_H
#define MARKER_H

#include <common.h>

struct s_marker: public s_card {
   Id Part;
   Id Point_mass;
   Vec3 Qp, Zp, Xp;
   Euler Reuler;
   Id Flex_body;
   Id Node_id;

   Boolean _Part, _Point_mass, _Qp, _Zp, _Xp, _Reuler;
   Boolean _Flex_body,_Node_id;
   Boolean _Floating;
   Boolean _UseXP;

   enum Deftype {
      ANGLES,
      POINT,
      EYEMATRIX
   };
   
   Deftype Mode;
   
   s_marker (void);
   ~s_marker (void);
   ostream& Print(ostream& out) const;
   void Translate(ostream& out);
   inline const char* const Gettype (void) const;
   Boolean Test();
   
   /* Only for marker section */
   void queue (void);
   void dequeue (void);
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch(param) {
	 case PART: Set_Param (Part,value,_Part); break;
	 case POINT_MASS : Set_Param (Point_mass, value, _Point_mass); break;
	 case QP : Set_Param (Qp, value, _Qp); break;
	 case XP : Set_Param (Xp, value, _Xp); break;
	 case ZP : Set_Param (Zp, value, _Zp); break;
	 case FLEX_BODY : Set_Param (Flex_body, value, _Flex_body) ; break;
	 case NODE_ID: Set_Param (Node_id, value, _Node_id); break;
	 case REULER: Set_Param (Reuler, value, _Reuler); break;
	 default : out_error(3,Find_Token(param)); break;
	}
	if (ES) out_error (1, Find_Token(param));
	return;
     }
};

#endif
