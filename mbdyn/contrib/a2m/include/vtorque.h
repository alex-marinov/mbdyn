#ifndef VTORQUE_H
#define VTORQUE_H

#include <common.h>

struct s_vtorque: public s_card {
   Id Node;
   Id Jfloat;
   Id Rm;
   double Tx,Ty,Tz;
   double Function;
   
   Boolean _Node,_Jfloat,_Rm,
           _Tx,_Ty,_Tz,_Function;
   
   s_vtorque(void);
   ~s_vtorque(void);
   inline const char* const Gettype (void) const;
   Boolean Test(void);
   ostream& Print (ostream& out) const;
   void Translate (ostream& out);
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch (param) {
	 case I: Set_Param (Node,value,_Node); break;
	 case RM: Set_Param (Rm,value,_Rm); break;
	 case JFLOAT: Set_Param (Jfloat,value,_Jfloat); break;
	 case TX: Set_Param (Tx,value,_Tx); break;
	 case TY: Set_Param (Ty,value,_Ty); break;
	 case TZ: Set_Param (Tz,value,_Tz); break;
	 case FUNCTION: Set_Param (Function,value,_Function); break;
	 default: out_error (3,Find_Token(param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
	return;
     }
};

#endif
