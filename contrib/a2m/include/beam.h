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

#ifndef BEAM_H
#define BEAM_H

#include <common.h>

struct s_beam: public s_card {
   Id Node[2];
   double Length;
   double Ixx,Iyy,Izz;
   double Area,Asy,Asz;
   double Gmodulus,Emodulus;
   double* Cmatrix;
   double Cratio;
   
   // Grandezze booleane di presenza dei parametri
   Boolean _Node1, _Node2, _Length, _Ixx, _Iyy, _Izz;
   Boolean _Area, _Asy, _Asz, _Gmodulus, _Emodulus, _Cmatrix, _Cratio;
   
   s_beam (void);
   ~s_beam (void);
   ostream& Print (ostream& out) const;
   void Translate(ostream& out);
   inline const char* const Gettype (void) const;
   Boolean Test();
   
   template <class T>
   void Set_Part (int param, T value)
   {	
      switch(param) {
       case I : Set_Param (Node[0], value, _Node1); break;
       case J : Set_Param (Node[1], value, _Node2); break;
       case LENGTH: Set_Param (Length, value, _Length); break;
       case IXX: Set_Param (Ixx, value, _Ixx); break;
       case IYY: Set_Param (Iyy, value, _Iyy); break;
       case IZZ: Set_Param (Izz, value, _Izz); break;
       case AREA: Set_Param (Area, value, _Area); break;
       case ASY : Set_Param (Asy, value, _Asy); break;
       case ASZ : Set_Param (Asz, value, _Asz); break;
       case EMODULUS: Set_Param (Emodulus, value, _Emodulus); break;
       case GMODULUS: Set_Param (Gmodulus, value, _Gmodulus); break;
       case CRATIO: Set_Param (Cratio, value, _Cratio); break;
       case CMATRIX: Set_Array (Cmatrix, value, 21, _Cmatrix); break;
       default: out_error(3,Find_Token(param)); break;
      }	
      if (ES) out_error (1, Find_Token(param));
      return;
   }
};

#endif
