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

       
	
