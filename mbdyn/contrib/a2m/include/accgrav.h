#ifndef ACCGRAV_H
#define ACCGRAV_H

#include <common.h>

struct s_accgrav : public s_card {
   double Igrav;
   double Jgrav;
   double Kgrav;
   
   Boolean _Igrav,_Jgrav,_Kgrav;
   
   s_accgrav(void);
   ~s_accgrav(void);
   ostream& Print (ostream& out) const;
   void Translate (ostream& out);
   inline const char* const Gettype (void) const;
   Boolean Test();
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch (param) {
	 case IGRAV: Set_Param (Igrav,value,_Igrav); break;
	 case JGRAV: Set_Param (Jgrav,value,_Jgrav); break;
	 case KGRAV: Set_Param (Kgrav,value,_Kgrav); break;
	 default : out_error (3, Find_Token(param)); break;
	}
	if (ES) out_error(1,Find_Token(param));
     }
};

#endif

