#ifndef VFORCE_H
#define VFORCE_H

#include <common.h>

struct s_vforce: public s_card {
   Id Node;
   Id Jfloat;
   Id Rm;
   double Fx,Fy,Fz;
   double Function;
   
   Boolean _Node,_Jfloat,_Rm,
           _Fx,_Fy,_Fz,_Function;
   
   s_vforce(void);
   ~s_vforce(void);
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
	 case FX: Set_Param (Fx,value,_Fx); break;
	 case FY: Set_Param (Fy,value,_Fy); break;
	 case FZ: Set_Param (Fz,value,_Fz); break;
	 case FUNCTION: Set_Param (Function,value,_Function); break;
	 default: out_error (3,Find_Token(param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
	return;
     }
};

#endif
