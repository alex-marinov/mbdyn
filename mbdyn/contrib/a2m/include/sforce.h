#ifndef FORCE_H
#define FORCE_H

#include <common.h>

struct s_sforce: public s_card {
   Id Node[2];
   Direction_Mode Mode;
   Boolean Actiononly;
   double Value;
   
   Boolean _Node1,_Node2,_Mode,
   _Value,_Actiononly;
   
   s_sforce(void);
   ~s_sforce(void);
   inline const char* const Gettype(void) const;
   Boolean Test(void);
   ostream& Print (ostream& out) const;
   void Translate (ostream& out);
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch (param) {
	 case I: Set_Param (Node[0],value,_Node1); break;
	 case J: Set_Param (Node[1],value,_Node2); break;
	 case _SET_:Set_Param(Mode,value,_Mode); break;
	 case ACTIONONLY:Set_Param(Actiononly, value,_Actiononly); break;
	 case FUNCTION: Set_Param(Value,value,_Value); break;
	 default: out_error (3,Find_Token(param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
	return;
     }
};

#endif

