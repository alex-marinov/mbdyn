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


#ifndef DEBUG_H
#define DEBUG_H

#include <storage.h>

// LIBRERIA CONTENENTE SEMPLICI FUNZIONI DI DEBUGGING

void DEBUGCOUT (const char* s);

template <class T>
void RESCUE ( T*& p, unsigned int N )
{
   char *commento = new char[160];
   sprintf (commento,"[DEBUG] : [%s OF LABEL %i] NOT DEFINED - RESCUE..",p->Gettype(),N); 
   DEBUGCOUT (commento);
   p = new T(N);
   return;
}

template <class T>
void DEBUG_AND_RESCUE (T*& p, unsigned int N )
{
   RESCUE (p,N);
   return;
}

template <class T>
void CHECK_AND_DEBUG (T* p, T*& p1, unsigned int N, MBDyn_deck& cd )
{
   
   /* Questa routine controlla se il puntatore alla card (p) è nullo
    * In questo caso crea una nuova card e assegna l'indirizzo a p1,
    * mentre se p è non nullo, rende uguali i due puntatori */
   
   char* s = new char[180];
   char* t = new char[180];
   if (p==NULL) {
     p1 = new T(N);
     sprintf (s,"\n# [DEBUG] - %s of label %i created by DEBUGGER",
	      p1->Gettype(),N);
     sprintf (t,"[DEBUG] : [UNREFERENCED CARD - CREATE NEW..]");
     p1->Remark(s);
     cout << t << endl;
     cd.insert(MBDyn_entry(N,(MBDyn_card*) p1));
   } else
     p1=p; /* la card esiste: rende uguali i due puntatori */
   return;
}
   
#endif
