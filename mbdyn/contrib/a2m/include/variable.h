#ifndef VARIABLE_H
#define VARIABLE_H

#include <common.h>

struct s_variable: public s_card {

   double Ic;
   char* Expression;
   double Value;
   
   Boolean _Ic,_Expression,_Value;
   
   s_variable(void);
   ~s_variable(void);
   inline const char* const Gettype(void) const;
   Boolean Test(void);
   ostream& Print (ostream& out) const;
   void Translate (ostream& out);
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch (param) {
	 case IC: Set_Param (Ic,value,_Ic); break;
	 case FUNCTION: Set_Param (Value,value,_Value); break;
	 default: out_error (3,Find_Token(param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
	return;
     }
};

#endif

