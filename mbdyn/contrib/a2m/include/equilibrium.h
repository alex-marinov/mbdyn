#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include <common.h>
#include <output.h>

struct s_equilibrium : public s_card {
   Angle Alimit;
   double Error;
   double Imbalance;
   Id Maxit;
   Boolean* Pattern;
   double Stability;
   double Tlimit;
   int ncoord_pattern;
   
   Boolean _Alimit, _Error, _Imbalance, _Maxit,
           _Pattern, _Stability, _Tlimit;
   
   s_equilibrium(void);
   ~s_equilibrium (void);
   ostream& Print(ostream& out) const;
   void Translate(ostream& out);
   inline const char* const Gettype (void) const;
   Boolean Test();
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch(param) {
	 case ALIMIT: Set_Param(Alimit,value,_Alimit); break;
	 case ERROR: Set_Param(Error, value,_Error); break;
	 case IMBALANCE : Set_Param(Imbalance,value,_Imbalance); break;
	 case MAXIT : Set_Param (Maxit,value,_Maxit); break;
	 case PATTERN: Set_Array (Pattern,value,ncoord_pattern,_Pattern); break;
	 case STABILITY: Set_Param (Stability,value,_Stability); break;
	 case TLIMIT : Set_Param (Tlimit,value,_Tlimit); break;
	}
	if (ES) out_error (1,Find_Token(param));
	return;
     }
};

#endif

   
