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

#ifndef STORAGE_H
#define STORAGE_H

#include <cards.h>

extern Boolean THROW_MODE;
extern Boolean OVERWRITE_MODE;
extern Boolean REMOVE_REMARK;

// DICHIARAZIONE DELLE STRUTTURE DATI PER IL CARD STORAGE

// STRUTTURE ADAMS
typedef map < Id,s_card*,less<Id> > card_deck;
typedef card_deck::value_type card_entry;
typedef card_deck::iterator p_card_entry;

// STRUTTURE TOKENS
typedef map < short int, short int, less<short int> > token_deck;
typedef token_deck::value_type token_entry;
typedef token_deck::iterator p_token_entry;

// STRUTTURE MBDYN
typedef map < Id, MBDyn_card*, less<Id> > MBDyn_deck;
typedef MBDyn_deck::value_type MBDyn_entry;
typedef MBDyn_deck::iterator p_MBDyn_entry;

// Struttura per il reference Marker To Part
/* la presente subroutine associa a ciascun marker
 * l'identificativo di una parte o di un pointmass */

struct MTP {
   Id Num;
   enum Type {
      PARTTYPE,
      POINTMASSTYPE,
      BEAMTYPE
   };
   Type El_type;
   MTP (Id N, MTP::Type q) : Num(N),El_type(q) {}
   MTP () : Num(0),El_type(PARTTYPE) {}
   ~MTP () {}
};

/* MTP associa a ogni Marker Adams una parte MBDYN */
typedef map < Id, MTP, less<Id> > MTP_deck;
typedef MTP_deck::value_type MTP_entry;
typedef MTP_deck::iterator p_MTP_entry;

/* Mappa per creare i riferimenti label to label */
typedef map < Id, Id, less<Id> > LTL_deck;
typedef LTL_deck::value_type LTL_entry;
typedef LTL_deck::iterator p_LTL_entry;

typedef LTL_deck PTR_deck,MTR_deck,PTN_deck,BTP_deck,BTR_deck;
typedef LTL_entry PTR_entry,MTR_entry,PTN_entry,BTP_entry,BTR_entry;
typedef p_LTL_entry p_PTR_entry,p_MTR_entry,p_PTN_entry,
  p_BTP_entry,p_BTR_entry;

/* STRUTTURE CHE SERVONO A GETFREE LABEL PER RESTITUIRE UN VALORE DI*/
/* LABEL NON ANCORA UTILIZZATO. I TRE FALDONI SONO ELEMENTI,NODI, E */
/* REFERENCE */

typedef map < Id, Boolean, less<Id> > Label_deck;
typedef Label_deck::value_type Label_entry;
typedef Label_deck::iterator p_Label_entry;

// FUNZIONE DI MEMORIZZAZIONE DELLE CARD DI ADAMS

template <class T>
void Set_Card (T*& pp, Id L, card_deck& cd)
{
   Boolean WRONG=N;
   Boolean TOINSERT=Y;
   s_card*& curr_card=(s_card*&) pp;
   p_card_entry idx;
   (pp)->label=L;
   idx=cd.find(L);
   WRONG=pp->Test(); 
   if ((*idx).second != NULL) {
    // CONTROLLA IL MODO DI SOVRASCRITTURA
      if (OVERWRITE_MODE==N) { 
        out_error (5,"");
        TOINSERT=N;
      }
   }
   if (WRONG==Y) {
     if (THROW_MODE==Y) TOINSERT=N;
   }
   if (TOINSERT==Y) cd.insert (card_entry(L, curr_card));
   pp = new T;
   return;
}


// FUNZIONE DI VISUALIZZAZIONE DEI CARD DECK DI ADAMS

template <class T>
void Display_deck (card_deck& cd, T*& pp, ostream& out)
{
   p_card_entry idx;
   for (idx=cd.begin(); idx != cd.end(); idx++) {
      pp = (T*) (*idx).second;
      pp->Print(out);
   }
   return;
}

// FUNZIONE DI VISUALIZZAZIONE DEI CARD DECK DI MBDYN

void Restart_MBDYNDeck (MBDyn_deck&, ostream&);

// FREE LABEL TO GET FROM A MBDYN DECK
Id GetFreeLabel (MBDyn_deck&);
Id GetFreeLabel (MBDyn_deck&,Id);
// FIND A CARD WITH LABEL ID IN A DECK
MBDyn_card* Find_MBCard (Id, MBDyn_deck&);
// FIND A CARD IN RECIPIENT
p_formula_entry trova (formula_map&,char*);


/* EXPERIMENTAL CODE BELOW*/
template <class T>
void Translate_deck (const char* msg, card_deck& cd, T* pp, ostream& out)
{
   T* temp;
   p_card_entry idx;
   if (cd.size()) { 
      cout << msg << endl;
      for (idx=cd.begin(); idx != cd.end(); idx++) {
	 temp = (T*) (*idx).second;
	 temp->Translate(out);
      }
   }
   return;
}

#endif
