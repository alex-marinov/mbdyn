
#ifndef SPRINGDAMPER_H
#define SPRINGDAMPER_H

#include <common.h>

struct s_springdamper: public s_card {
   Id Node1;
   Id Node2;
   Direction_Mode Mode;
   double C,K,Force,Length;
   double Ct,Kt,Torque;
   Angle ANgle;
   
   Boolean _Node1,_Node2, _Mode,
           _C,_K,_Force,_Length,
           _Ct,_Kt,_Torque,_ANgle;
   
   s_springdamper(void);
   ~s_springdamper (void);
   inline const char* const Gettype(void) const;
   Boolean Test(void);
   ostream& Print (ostream& out) const;
   void Translate (ostream& out);
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch(param) {
	 case I:Set_Param(Node1,value,_Node1); break;
	 case J:Set_Param(Node2,value,_Node2); break;
	 case C0: Set_Param (C, value, _C); break;
	 case K0: Set_Param (K, value, _K); break;
	 case FORCE: Set_Param (Force, value, _Force); break;
	 case LENGTH: Set_Param (Length, value, _Length); break;
	 case CT: Set_Param (Ct, value, _Ct); break;
	 case KT: Set_Param (Kt, value, _Kt); break;
	 case TORQUE: Set_Param (Torque, value, _Torque); break;
	 case ANGLE: Set_Param (ANgle, value, _ANgle); break;
	 case _SET_: Set_Param (Mode, value, _Mode); break;
	 default: out_error(3,Find_Token (param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
	return;
     }
};

#endif
