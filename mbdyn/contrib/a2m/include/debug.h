
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
