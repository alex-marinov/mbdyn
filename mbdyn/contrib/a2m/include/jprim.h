#ifndef JPRIM_H
#define JPRIM_H

#include <common.h>

struct s_jprim: public s_card {
   Id Node1;
   Id Node2;
   Joint_Primitive Type;
   
   Boolean _Node1,_Node2,_Type;
   
   s_jprim(void);
   ~s_jprim(void);
   ostream& Print (ostream& out) const;
   void Translate(ostream& out);
   inline const char* const Gettype (void) const;
   Boolean Test();
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch (param) {
	 case I: Set_Param (Node1, value, _Node1); break;
	 case J: Set_Param (Node2, value, _Node2); break;
	 case _SET_ : Set_Param (Type, value, _Type); break;
	 default: out_error (3,Find_Token(param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
	return;
     }
};

#endif

