#ifndef CARDS_H
#define CARDS_H

// ROUTINE DI FONDAMENTALE GESTIONE DELLE CARDS
#include <defs.h>
#include <errrec.h>
#include <mbdyn.h>
#include <mathem.h>
#include <iostream.h>
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
