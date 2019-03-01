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

//                               STORAGE.CC                                   

#include <storage.h>

/* GETFREELABEL: Restituisce un valore di Label non ancora utilizzato */

Id GetFreeLabel (MBDyn_deck& L)
{
   Id Ref1,Ref2;
   Id ValRet;
   Boolean STATE=N;
   Boolean EMPTY=N;
   int i=1;
   p_MBDyn_entry PL1,PL2;
   if (L.empty()) EMPTY=Y; 
   if (!EMPTY) {
      PL1=L.begin();
      while (STATE==N) {
	 if (PL2==L.end()) {
	    STATE=Y;
	    ValRet=(Ref2+1);
	 }
	 Ref1=(*PL1).first;
	 PL2=PL1;
	 Ref2=(*(++PL2)).first;
	 PL1++;
	 if ((Ref2-Ref1) > 1) {
	    STATE=Y;
	    ValRet=(Ref1+1);
	 }
      }
   }
   else {
      /* Essendo la lista vuota assegna come valore il primo utile,cioè 1*/
      ValRet=1;
   }
   return ValRet;
}

Id GetFreeLabel (MBDyn_deck& TANK, Id PREF)
{
   //
   /* questa funzione restituisce un valore di label libero per il
    * contenitore TANK. Se il valore PREF non risulta stato utilizzato
    * ancora, utilizza PREF, altrimenti ritorna una chiamata a GetFreeLabel
    * per una label non ancora utilizzata */
   //
   MBDyn_card* TEST;
   Id TORETURN;
   TEST=Find_MBCard (PREF,TANK);
   if (TEST==NULL) TORETURN=PREF;
   else TORETURN= GetFreeLabel (TANK);
   return TORETURN;
}

MBDyn_card* Find_MBCard (Id L, MBDyn_deck& TANK)
{
   Id Check;
   MBDyn_card* Object;
   p_MBDyn_entry Iterator;
   Iterator=TANK.find (L);
   Object = (*Iterator).second;
   Check=(*Iterator).first;
   if (Check!=L) return NULL;
   else return Object;
}  

void Restart_MBDYNDeck (MBDyn_deck& cd, ostream& out)
{
   p_MBDyn_entry idx;
   MBDyn_card* CARD;
   for (idx=cd.begin(); idx != cd.end(); idx++) {
      CARD = (*idx).second;
      if ((CARD->_remark_ != NULL) & REMOVE_REMARK==N) 
	out << CARD->_remark_ << endl;
      CARD->Restart(out);
   }
   return;
}

p_formula_entry trova(formula_map& m,char* key)
{
   p_formula_entry idx;
   for (idx=m.begin();idx!=m.end();idx++)
     if (strcmp( (*idx).first,key ) == 0 ) {
	return idx;
     }
   idx=m.find(key);
   return idx;
}

       
	
