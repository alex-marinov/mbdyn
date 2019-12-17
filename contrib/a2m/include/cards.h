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


#ifndef CARDS_H
#define CARDS_H

// ROUTINE DI FONDAMENTALE GESTIONE DELLE CARDS
#include <defs.h>
#include <errrec.h>
#include <mbdyn_a2m.h>
#include <mathem.h>
#include <iostream>
#include <stdio.h>
#include <output.h>

// GENERIC CARD DECLARATION

struct s_card {
   Id label;
   char* _remark_;
   enum {
      _BEAM,
	_MARKER,
	_PART,
        _JOINT,
      
	_LAST_CARD
   };
   s_card(void);
   virtual ostream& Print (ostream& out) const=0;
   virtual void Translate (ostream& out)=0;
   virtual inline const char* const Gettype(void) const=0;
   virtual Boolean Test ()=0;

   void Comment (char*);
   void Display_Formula(ostream&) const;
   void Store_Formula(char*,char*);
   formula_map Recipient;
};

// ARRAY AND PARAMETERS STORING PROCEDURES FOR ALL CARDS

template <class T>
void Set_Array (T& param, T value, int nc, Boolean& Flag)
{
   ES=N;
   if (Flag==N) {
      for (int k=0; k<nc; k++) param[k]=value[k];
      Flag=Y;
      return;
   }
   ES=Y;
   return;
}

template <class T,class Q>
void Set_Array (T& param, Q value, int nc, Boolean& Flag)
{
   out_error (4,"ARRAY TYPE");
   return;
}

template <class T>
void Set_Param (T& param, T value, Boolean& Flag)
{
   ES=N;
   if (Flag==N) {
      param=value;
      Flag=Y;
      return;
   } 
   ES=Y;
   return;
}

template <class T,class Q>
void Set_Param (Q& param, T value, Boolean& Flag)
{
   out_error(4,"PARAM TYPE");
   return;
}

#endif
